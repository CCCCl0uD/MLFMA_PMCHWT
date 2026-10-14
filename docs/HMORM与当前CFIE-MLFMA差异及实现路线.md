# HMORM 与当前 CFIE-MLFMA 差异及实现路线

本文档说明论文 HMORM 与当前工程中 CFIE-MLFMA 的主要区别，并拆解“物理层降阶”和“算法层加速”如何实现。重点结合当前代码结构，说明如果要在 `PMCHWT_MLFMA/MLFMA` 工程中落地 HMORM，应改哪些数据结构、预处理流程、矩阵向量乘法和 RCS 输出流程。

参考文档：

- `docs/HMORM公式原理拆解与应用扩展.md`
- 源论文：Si Jia Guo, Yu Mao Wu, Yuan An Liu, "An Efficient Hierarchical Model Order Reduction Method for Analyzing Electromagnetic Scattering from Large-Scale Periodic Arrays", IEEE TAP, DOI: 10.1109/TAP.2026.3720418

## 1. 一句话结论

当前代码中的 CFIE-MLFMA 是：

```text
CFIE 物理方程 + RWG 未知量 + RWG 级 MLFMA 矩阵向量乘法
```

论文中的 HMORM 是：

```text
CFIE 物理方程 + CBF 降阶未知量 + CBF 级 MLFMM 矩阵向量乘法
```

所以二者最大的区别不是 CFIE 公式本身，而是求解空间和 MLFMA 叶层对象不同：

| 项目 | 当前 CFIE-MLFMA | HMORM |
|---|---|---|
| 积分方程 | CFIE | 仍可用 CFIE |
| 基函数 | RWG basis | CBF，即 RWG 的局部线性组合 |
| 未知量 | 每条 RWG 边一个电流系数 `I_i` | 每个块少量 CBF 系数 `alpha_{b,k}` |
| 求解器维度 | `N_RWG` | `K_total = sum_b K_b` |
| 叶层聚合 | 对每条 RWG 做 `Vsmi[id] * x[id]` | 对每个 CBF 做 `S_{b,k} * alpha_{b,k}` |
| 叶层解聚合/测试 | 对每条 RWG 用 `Vfmj[rwgid]` 测试 | 对每个 CBF 用 `V_{b,j}` 测试 |
| 近场矩阵 | RWG-RWG 的 `Z_near` | CBF-CBF 的 `Z_near_red` |
| RHS | RWG 级 `Vm` | 投影后的 CBF 级 `Vm_red = Phi^T Vm` |
| 输出电流 | 直接得到 RWG 电流 | 先解 CBF 系数，再回投 `I_RWG = Phi alpha` |
| 加速来源 | 只加速远场耦合 | 同时减少未知量、近场规模、远场聚合规模和 Krylov 维度 |

## 2. 当前代码中的 CFIE-MLFMA 是怎样工作的

### 2.1 构造阶段

当前 `MLFMM` 构造函数在 `Code/MLFMM.cpp` 中完成预处理。对 PEC CFIE 分支，主要步骤是：

1. 调用 `AlgoMLFMM::computeBase_MLFMM(...)` 计算层级角谱、插值阶数、`kp` 网格、积分权重，并设置 `row = rwgs.size()`。
2. 调用 `AlgoMLFMM::computeTransferFactorMatrix(...)` 计算每层远场转移因子 `TF`。
3. 调用 `AlgoMLFMM::computeInterpolationMatrix(...)` 计算层间插值矩阵。
4. 调用 `ProcessFMM::computeVsmi_Vfmj_Cfie(...)` 计算每条 RWG 的源侧聚合函数 `Vsmi` 和测试侧函数 `Vfmj`。
5. 调用 `Zmartix::computeZ_fmm_pec_cfie(...)` 计算最细层近邻盒子之间的 RWG-RWG 近场矩阵 `Z_near`。
6. 根据单站或双站 RCS 配置进入 `mlfmm_Mono_Pec_Cfie` 或 `mlfmm_Dual_Pec_Cfie`。

核心入口可对应到这些代码位置：

- `Code/MLFMM.cpp`: 构造函数中 `computeBase_MLFMM`、`computeTransferFactorMatrix`、`computeVsmi_Vfmj_Cfie`、`computeZ_fmm_pec_cfie`。
- `Code/AlgoMLFMM.h`: `computeBase_MLFMM` 中直接令 `row = rwgSize`，说明当前全局未知量数就是 RWG 数。
- `Code/ProcessFMM.h`: `computeVsmi_Vfmj_Cfie` 为每条 RWG 生成 MLFMA 叶层谱函数。
- `Code/Zmartix.h`: `computeZ_fmm_pec_cfie` 为每条 RWG 测试函数存储邻近 RWG 源基函数的阻抗条目。

### 2.2 求解阶段

`MLFMM::matrix_solver(...)` 中定义了矩阵向量乘法 `ax`。PEC CFIE 情况下，`x` 的长度为 `row`，也就是全局 RWG 数。每次 GMRES 迭代都会执行：

```text
x_RWG
  -> mexp: 叶层 RWG 源聚合
  -> m2m: 逐层向上聚合
  -> m2l: 远场组间转移
  -> l2l: 逐层向下解聚合
  -> leaf testing: 每条 RWG 用 Vfmj 做测试
  -> add near field: 加 RWG-RWG Z_near
  -> w_RWG = A_RWG x_RWG
```

当前代码中最关键的片段是：

```cpp
ProcessMLFMM::clearSmBm(Sm, Bm);
ProcessMLFMM::mexp(Sm, octreeNodes_[maxlvl], Vsmi, kp_lvl_k1[0], x, 0, maxlvl);
ProcessMLFMM::m2m(Sm, octreeNodes_, kp_lvl_k1, interpol_k1, maxlvl, levelSpan, wave.k1());
ProcessMLFMM::m2l(Sm, Bm, octreeNodes_, TF, kp_lvl_k1, WGL_k1, WGL_phi_k1, phi_level_k1, maxlvl, levelSpan);
ProcessMLFMM::l2l(Bm, octreeNodes_, kp_lvl_k1, phi_level_k1, interpol_k1, maxlvl, levelSpan, wave.k1());
```

然后在最细层：

```cpp
for (int r = 0; r < rwgNum; ++r) {
    int rwgid = node->rwgIndices[r]->rwgid;
    const auto& V = Vfmj[rwgid];
    ...
    val += ProcessMLFMM::lexpKpDot(B[tt * pN_k1 + pp], V[tt * pN_k1 + pp]);
    ...
    val += x[Z_near_id[rwgid][j]] * Z_near[rwgid][j];
    w[rwgid] = val;
}
```

这说明当前代码的远场和近场最终都回到 RWG 行上，求解器看到的是 `A_RWG x_RWG = b_RWG`。

### 2.3 RHS 与 RCS

CFIE 激励向量在 `Code/Vm.h` 中由

```cpp
Vm[i] = alpha_ * V1[i] + (1.0 - alpha_) * V2[i];
```

得到，其中 `V1` 是 EFIE 右端项，`V2` 是 MFIE 右端项。

单站 RCS 路径在 `Code/RCS.h` 的 `computeMonoStatic(...)` 中逐角度执行：

```text
computeV(solver, kinc, einc, hinc)
solver.matrix_solver(solver.row, I_solve, Vec_R, ...)
RCSUtils::calculateRCS(solver, I_solve, ...)
```

所以当前实现每个角度解的是 RWG 电流 `I_solve[0:row)`，RCS 计算也直接遍历所有 RWG。

## 3. HMORM 的“物理层降阶”如何实现

物理层降阶的目标是把每个局部块内很多 RWG 自由度压缩为少量能代表主要电流形态的 CBF。

### 3.1 数学形式

当前 CFIE-MLFMA 的未知量是

```math
J(r') = \sum_{i=1}^{N_{\mathrm{RWG}}} I_i f_i(r')
```

HMORM 改为按块展开：

```math
J(r') =
\sum_{b=1}^{M}\sum_{k=1}^{K_b}
\alpha_{b,k} g_{b,k}(r')
```

其中每个 CBF 又是本块 RWG 的线性组合：

```math
g_{b,k}(r') =
\sum_{n=1}^{N_b}
\Phi_b(n,k) f_{b,n}(r')
```

这里 `Phi_b` 就是代码中需要保存的 CBF 系数矩阵，维度为：

```text
N_RWG_in_block x K_block
```

于是全局 RWG 电流和 CBF 系数的关系是：

```math
I_{\mathrm{RWG}} = \Phi \alpha_{\mathrm{CBF}}
```

其中 `Phi` 是按块对角拼起来的稀疏/块稠密矩阵。

### 3.2 CBF 怎么生成

论文采用的是 SVD-based MPW CBF：

1. 选一个原始块，例如一个周期阵元。
2. 如果阵元电连接或块边界处电流连续性重要，则构造扩展块。
3. 对扩展块施加多组平面波激励，得到多个 RHS。
4. 解局部 CFIE 方程：

```math
Z_{\mathrm{Ext},b}^{SI}
I_{\mathrm{Ext},b}^{MPW}
=
V_{\mathrm{Ext},b}^{MPW}
```

5. 对响应矩阵做 SVD：

```math
I_{\mathrm{Ext},b}^{MPW}
=
U_b S_b V_b^H
```

6. 按阈值 `eta0` 保留主导左奇异向量，得到扩展块 CBF。
7. 用选择矩阵 `C_b` 从扩展块中截回原始块：

```math
\Phi_b =
I_{\mathrm{Ori},b}^{CBF}
=
C_b I_{\mathrm{Ext},b}^{CBF}
```

### 3.3 CBF 中的 RWG 是如何选取的

CBF 中的 RWG 不是重新生成的一套网格，也不是从全局 RWG 中随机挑选若干条边。它的选取遵循“先按物理块选 RWG，再按扩展块补邻域 RWG，最后通过 SVD 形成这些 RWG 的线性组合”这一逻辑。

对第 `b` 个 CBF block，需要区分两组 RWG：

| 名称 | 记号 | 作用 | 是否进入最终全局未知量 |
|---|---|---|---|
| 原始块 RWG | `J_{\mathrm{Ori}}^b` 或 `rwgIds` | 属于第 `b` 个阵元/单元/leaf block 的 RWG | 是 |
| 扩展块 RWG | `J_{\mathrm{Ext}}^b` 或 `extRwgIds` | 原始块加上邻域内用于捕捉边界连续电流和近耦合的 RWG | 不全部进入，最后只截回原始块 |

二者满足：

```math
J_{\mathrm{Ori}}^b \subset J_{\mathrm{Ext}}^b
```

#### 3.3.1 原始块 RWG 的选取

原始块 RWG 是最终 CBF 所覆盖的自由度。对周期阵列，最自然的定义是：

```text
一个周期单元内的所有 RWG = 一个原始块的 RWG
```

如果一个 RWG 的两个三角形都位于同一单元内，它显然属于该单元。麻烦在于边界 RWG，尤其是电连接阵列中跨越相邻单元边界的 RWG。推荐规则是：

1. 用 RWG 边中心或两个三角形中心的几何位置判断所属单元。
2. 若 RWG 跨单元边界，则按边中心所在单元归属到一个原始块，不能在两个原始块中重复计数。
3. 所有全局 RWG 必须且只能属于一个原始块，否则回投 `I_{\mathrm{RWG}} = \Phi\alpha` 时会出现漏计或重复叠加。

也就是说，原始块 RWG 是一个全局剖分：

```math
\bigcup_{b=1}^{M} J_{\mathrm{Ori}}^b
=
\{1,2,\ldots,N_{\mathrm{RWG}}\},
\qquad
J_{\mathrm{Ori}}^a \cap J_{\mathrm{Ori}}^b = \varnothing,\ a\ne b
```

对当前代码，如果先做快速验证版，可以直接令：

```text
J_Ori^b = octreeNodes_[maxLevel_][b]->rwgIndices
```

也就是把当前 MLFMA 最细层 leaf node 中的 RWG 当作一个 CBF block 的原始 RWG。这样最容易跑通，但它未必等价于论文中的“一个周期单元一个 CBF block”。

正式周期阵列版应优先按阵元或周期单元编号来选 RWG，而不是按 octree 自动分桶来选。理想结果是：

```text
unitCellId -> rwgIds -> CBFBlock
CBFBlock  -> 对应一个 finest-level cube
```

#### 3.3.2 扩展块 RWG 的选取

扩展块 RWG 用于生成 CBF 时的局部局域问题。它的目标不是增加最终未知量，而是让局部响应矩阵能“看见”邻域耦合。

对第 `b` 个原始块，扩展块 RWG 可按以下规则选：

```text
extRwgIds(b) =
    原始块 rwgIds(b)
  + 距离原始块几何边界不超过 lExt 的邻域 RWG
```

其中 `lExt` 是扩展长度。论文给出的范围是 `0.025 lambda` 到 `0.4 lambda`。对电连接贴片阵列，`0.025 lambda` 可理解为向外扩展一圈相邻三角面片。

在实现上有三种常见判据：

1. 按单元邻接扩展：把上下左右或 8 邻域单元中的 RWG 加入扩展块。
2. 按几何距离扩展：如果 RWG 边中心、三角形中心或三角形任一点距离原始块包围盒小于 `lExt`，则加入。
3. 按网格拓扑扩展：从原始块边界 RWG 出发，沿共享顶点或共享边扩展一圈或多圈三角面片。

电连接阵列优先推荐第 3 种或“第 2 种 + 拓扑修正”，因为电流连续性沿网格连通关系传播；单纯按包围盒距离可能会把几何上近但电流不连通的面片纳入。

#### 3.3.3 SVD 中每一行对应一个扩展块 RWG

局部响应矩阵

```math
I_{\mathrm{Ext},b}^{MPW}
```

的行对应扩展块中的 RWG，列对应不同平面波激励：

```text
rows    = extRwgIds.size()
columns = N_MPW
```

第 `i` 行就是 `extRwgIds[i]` 这条 RWG 在所有 MPW 激励下的电流响应。SVD 后保留的左奇异向量同样定义在扩展块 RWG 行空间上。

因此，一个 CBF 模态不是“选中了几条 RWG”，而是“在扩展块所有 RWG 上给出一组系数”：

```math
g_{b,k}^{Ext}(r')
=
\sum_{i=1}^{N_{\mathrm{Ext}}^b}
U_b(i,k) f_{\mathrm{extRwgIds}[i]}(r')
```

然后再用选择矩阵 `C_b` 截回原始块：

```math
g_{b,k}^{Ori}(r')
=
\sum_{i:\ \mathrm{extRwgIds}[i]\in J_{\mathrm{Ori}}^b}
U_b(i,k) f_{\mathrm{extRwgIds}[i]}(r')
```

最终保存到 `Phi_b` 的，只是原始块 RWG 对应的那些行：

```text
Phi_b(localOriRwg, k)
```

#### 3.3.4 哪些 RWG 不应该放进最终 CBF

以下 RWG 不应进入某个块的最终 `Phi_b`：

- 只属于邻居扩展区域、但不属于该原始块的 RWG；
- 已经归属到其他原始块的跨边界 RWG；
- 距离上接近但拓扑或材料区域不相关的 RWG；
- 在 PMCHWT 中属于另一种未知量却没有同步进入联合基的 RWG。

这些 RWG 可以参与扩展块 CBF 生成，但最终必须通过 `C_b` 截掉。否则一个物理 RWG 会被多个 CBF block 同时表示，导致全局自由度重复。

#### 3.3.5 对多个相同目标 RCS 的 RWG 选取

如果把 HMORM 用于多个相同目标 RCS，原始块 RWG 的定义通常更简单：

```text
一个目标实例上的所有 RWG = 一个原始块的 RWG
```

若单个目标电尺寸较大，应先把单目标再切成多个子块：

```text
一个目标实例 = 若干 CBFBlock
多个目标实例复用同一套 block 类型
```

扩展块可按目标间距离决定：

- 目标相距较远：`extRwgIds = rwgIds` 即可；
- 目标靠得很近、有强近场耦合：把相邻目标表面靠近区域的 RWG 加入扩展块；
- 目标接触或电连接：必须按拓扑连通或一圈面片扩展，不能只按目标编号硬切。

### 3.4 物理意义

CBF 不是凭空构造的数学基，而是“这个局部块在多方向、多极化入射下可能出现的主要电流模式”。因此它带有物理信息：

- 低阶 CBF 对应最容易被激发的主电流分布；
- 高阶奇异向量通常对应较弱、较局部或数值噪声占比较高的模式；
- SVD 截断相当于舍弃对整体散射贡献较小的局部电流自由度；
- 扩展块使 CBF 生成时见到邻域结构，能更好表示电连接边界处的连续电流。

这就是“物理层降阶”：不是简单删除网格，也不是粗网格化，而是在保持原始 RWG 网格的前提下，把解限制在一个物理响应子空间中。

### 3.5 相较当前 CFIE-MLFMA 的改进

当前 CFIE-MLFMA 没有物理层降阶：

```text
每条 RWG 边 = 一个未知量
```

HMORM 变成：

```text
每个块 = K 个 CBF 未知量，K 通常远小于该块 RWG 数
```

因此改进是全链条的：

- GMRES 向量长度变小；
- Krylov 正交化成本下降；
- 近场矩阵规模下降；
- 远场叶层聚合对象减少；
- 多 RHS 或多角度扫描时，CBF 和预处理数据可复用。

## 4. HMORM 的“算法层加速”如何实现

算法层加速不是重新发明 MLFMA，而是把现有 MLFMA 从 RWG 空间改成 CBF 空间。

### 4.1 当前代码的叶层聚合

当前 `ProcessMLFMM::mexp(...)` 做的是：

```cpp
Sm[leaf][node][ki] += Vsmi[rwgid][ki] * x[rwgid];
```

数学上对应：

```math
S_m(\hat{k}) =
\sum_{q\in G_m}
F_q(\hat{k}) I_q^{RWG}
```

也就是一个叶节点内所有 RWG 源基函数的远场谱叠加。

### 4.2 HMORM 的叶层聚合

HMORM 先预计算每个 CBF 的源谱：

```math
S_{b,k}(\hat{k})
=
\sum_{q=1}^{N_b}
\Phi_b(q,k) F_q(\hat{k})
```

然后每次迭代只做：

```math
S_b^{iter}(\hat{k})
=
\sum_{k=1}^{K_b}
S_{b,k}(\hat{k}) \alpha_{b,k}
```

和当前代码相比，区别是：

```text
当前：每次迭代遍历 node->rwgIndices
HMORM：每次迭代遍历 block->cbfModes
```

如果一个叶块有 `N_b = 300` 条 RWG，压缩后 `K_b = 30` 个 CBF，则叶层聚合循环理论上减少约 10 倍。

### 4.3 当前代码的叶层测试

当前代码在最细层对每条 RWG 做：

```cpp
val += dot(Bm_leaf[node][kp], Vfmj[rwgid][kp]);
```

数学上对应：

```math
w_p^{far}
=
\int R_p(\hat{k})\cdot B_m(\hat{k})d^2\hat{k}
```

结果是 RWG 行 `w[p]`。

### 4.4 HMORM 的叶层测试

HMORM 预计算每个 CBF 的测试函数：

```math
V_{b,j}(\hat{k})
=
\sum_{p=1}^{N_b}
\Phi_b(p,j) R_p(\hat{k})
```

远场测试直接得到 CBF 行：

```math
w_{b,j}^{far}
=
\int V_{b,j}(\hat{k})\cdot B_b(\hat{k})d^2\hat{k}
```

也就是矩阵向量乘法的输出从

```text
w_RWG[global_rwg_id]
```

变成

```text
w_CBF[global_cbf_id]
```

### 4.5 当前代码的近场项

当前 `Zmartix::computeZ_fmm_pec_cfie(...)` 为每条 RWG 存储其所在叶盒和近邻叶盒中的 RWG-RWG 阻抗：

```cpp
Z_near[rwgfid][count] =
    alpha * Zmn_near_E
  + (1.0 - alpha) * wave.eta1() * Zmn_near_M;
Z_near_id[rwgfid][count] = rwgSid;
```

每次矩阵向量乘时做：

```cpp
val += x[Z_near_id[rwgid][j]] * Z_near[rwgid][j];
```

这是 RWG 级近场。

### 4.6 HMORM 的近场项

HMORM 近场矩阵应先投影：

```math
Z_{ab}^{CBF}
=
\Phi_a^T Z_{ab}^{RWG} \Phi_b
```

其中 `a` 是测试块，`b` 是源块。矩阵尺寸从

```text
N_a x N_b
```

变成

```text
K_a x K_b
```

迭代中只做：

```math
w_a^{near}
=
\sum_{b\in Near(a)}
Z_{ab}^{CBF}\alpha_b
```

如果周期阵列严格重复，还可以按相对位移复用 `Z_{ab}^{CBF}`，不必每个块对都存一份。

### 4.7 层级传播仍可复用

`m2m`、`m2l`、`l2l` 的核心几何传播逻辑可以基本复用，因为它们只关心：

- 每层节点；
- 节点中心；
- 插值矩阵；
- 转移因子；
- `Sm/Bm` 的角谱数组。

真正需要变化的是最细层的“源谱如何生成”和“场谱如何测试”。这对当前代码是个好消息：不用整体推翻 MLFMA，只要把 leaf-level 的输入输出从 RWG 改成 CBF。

## 5. 为什么 HMORM 比当前 CFIE-MLFMA 更快、更省内存

设：

- `M` = 块数；
- `N` = 每块 RWG 数；
- `K` = 每块 CBF 数；
- `K << N`。

当前 CFIE-MLFMA 的主要迭代复杂度约为：

```math
O(MN\log(MN))
```

HMORM 改为：

```math
O(MK\log(MK))
```

这不是只省了一处，而是几处同时下降：

| 模块 | 当前 CFIE-MLFMA | HMORM 改进 |
|---|---|---|
| GMRES 向量 | `MN` | `MK` |
| Arnoldi 正交化 | 与 `MN` 成正比 | 与 `MK` 成正比 |
| 叶层聚合 | 每次迭代遍历 RWG | 每次迭代遍历 CBF |
| 叶层测试 | 输出 RWG 行 | 输出 CBF 行 |
| 近场矩阵 | RWG-RWG 条目 | CBF-CBF 小矩阵 |
| 近场存储 | 邻域内所有 RWG 对 | 邻域块对的降阶矩阵 |
| 多角度 RHS | 每角度 RWG RHS 和 RWG 解 | CBF/RWG 投影数据可复用 |

论文中大规模阵列能比 MLFMM 快十几倍，原因就是普通 MLFMM 只加速远场，但未知量仍然是 RWG；HMORM 把“参与远场加速的最小对象”也变成了 CBF。

## 6. 当前代码要实现 HMORM 的核心改造点

### 6.1 新增 CBF 数据结构

建议新增 `Code/HMORM.h` 和 `Code/HMORM.cpp`，第一版先只做 PEC CFIE。

核心结构可设计为：

```cpp
struct CBFBlock {
    int blockId = -1;
    int nodeId = -1;
    std::vector<int> rwgIds;       // 原始块 RWG
    std::vector<int> extRwgIds;    // 扩展块 RWG，可第一版为空或等于 rwgIds
    int K = 0;

    // Phi[col-major or row-major]:
    // size = rwgIds.size() * K
    // Phi(localRwg, mode)
    std::vector<std::complex<double>> Phi;

    // CBF source aggregation and testing functions.
    // size = K x kpCount
    std::vector<std::vector<Complex3D>> S_cbf;
    std::vector<std::vector<Complex3D>> V_cbf;
};

struct HMORMContext {
    std::vector<CBFBlock> blocks;
    std::vector<int> cbfOffset; // global alpha offset per block
    int totalCbf = 0;
};
```

第一版可让 `blockId == finest node index`，即每个 leaf node 是一个 CBF block。后续如果要严格复现论文，需要让 leaf cube 与周期单元一一对应，而不是完全依赖当前 octree 的自动边中心分桶。

### 6.2 叶节点与物理块同构

论文 HMORM 的关键是：

```text
finest-level cube <-> CBF block <-> 一个周期单元
```

当前 `OCTree::buildLeafNodes(...)` 是按 RWG 边中心落入 leaf cube：

```cpp
const auto& p = rwg.edgeCenterCenter;
int ix = static_cast<int>((p.x - globalBBoxMin_.x) / leafLen_);
...
nodePtr->rwgIndices.push_back(...)
```

这未必保证一个周期单元正好一个 leaf node。要实现论文式 HMORM，有两种选择：

1. 快速验证版：沿用当前 leaf node，当作 CBF block。优点是改动小；缺点是无法充分利用周期单元复用。
2. 正式周期阵列版：新增基于阵元编号或单元坐标的分块逻辑，使每个单元的 RWG 明确归入同一个 `CBFBlock`，再反过来约束或生成 MLFMM leaf node。

推荐路线是先做 1，再做 2。先让降阶矩阵向量乘跑通，比一开始就重构 octree 更稳。

### 6.3 生成 CBF

第一版建议先实现“非扩展块 CBF”，也就是 `extRwgIds = rwgIds`。跑通后再做扩展块。

需要实现：

```text
for each block b:
    build local dense CFIE Zbb using existing EFIE/MFIE operators
    build MPW RHS matrix Vmpw
    solve Zbb * I_MPWs = Vmpw
    SVD(I_MPWs)
    choose K by eta0
    Phi_b = leading left singular vectors
```

其中 `rwgIds` 和 `extRwgIds` 的构造建议先单独封装，避免后面投影 RHS、近场矩阵、远场谱函数时重复写映射逻辑：

```cpp
struct BlockRwgSelection {
    std::vector<int> rwgIds;       // 原始块 RWG，全局唯一归属
    std::vector<int> extRwgIds;    // 扩展块 RWG，可包含邻居 RWG
    std::vector<int> oriRowsInExt; // 原始块 RWG 在 extRwgIds 中的行号
};
```

快速验证版可以从当前 leaf node 直接生成：

```cpp
BlockRwgSelection selectRwgsFromLeafNode(const OCTree::Node* node)
{
    BlockRwgSelection sel;
    for (const RWGBase* rwg : node->rwgIndices) {
        sel.rwgIds.push_back(rwg->rwgid);
        sel.extRwgIds.push_back(rwg->rwgid);
        sel.oriRowsInExt.push_back((int)sel.oriRowsInExt.size());
    }
    return sel;
}
```

扩展块版本需要先建立 `rwgid -> originalBlockId` 和 `rwgid -> localOriIndex`，然后按几何距离或网格拓扑加入邻域 RWG：

```text
for each original block b:
    rwgIds = all RWGs whose owner block is b
    extRwgIds = rwgIds
    for each candidate RWG near block b:
        if distance/topology criterion <= lExt:
            append candidate RWG to extRwgIds if not duplicated
    oriRowsInExt[i] = row index of rwgIds[i] inside extRwgIds
```

`oriRowsInExt` 就是论文选择矩阵 `C_b` 的稀疏实现。实际代码里不需要显式存一个 0/1 矩阵，只要用这些行号从 `U_b` 中抽取：

```cpp
for (int i = 0; i < (int)sel.rwgIds.size(); ++i) {
    int rowInExt = sel.oriRowsInExt[i];
    for (int k = 0; k < K; ++k) {
        Phi[i * K + k] = U[rowInExt * Kfull + k];
    }
}
```

可复用当前代码：

- 局部 `Zbb` 可复用 `EFIE::computeEFIE_Zij` 和 `MFIE::computeMFIE_Zij` 的组合方式。
- 局部 RHS 可参考 `RHS::computeV_E`、`RHS::computeV_M`、`RHS::computeV_C`，但要支持多平面波方向和只对 block 内 RWG 计算。
- SVD 需要引入线性代数库。若当前工程未接 Eigen/LAPACK，可先用 Eigen 的 `JacobiSVD` 做小块验证。

CFIE 局部矩阵仍是：

```math
Z_{ij}^{CFIE}
=
\alpha Z_{ij}^{EFIE}
+
(1-\alpha)\eta Z_{ij}^{MFIE}
```

这和当前 `Zmartix::computeZ_fmm_pec_cfie(...)` 中的组合方式一致。

### 6.4 投影 RHS

当前 RHS 是 `Vm[rwgid]`。HMORM 需要生成 CBF RHS：

```math
V_b^{CBF}
=
\Phi_b^T V_b^{RWG}
```

代码上可新增：

```cpp
void projectRhsToCBF(
    const HMORMContext& ctx,
    const std::vector<std::complex<double>>& rhsRwg,
    std::vector<std::complex<double>>& rhsCbf)
{
    rhsCbf.assign(ctx.totalCbf, 0.0);
    for (const auto& b : ctx.blocks) {
        int off = ctx.cbfOffset[b.blockId];
        for (int k = 0; k < b.K; ++k) {
            std::complex<double> sum = 0.0;
            for (int i = 0; i < (int)b.rwgIds.size(); ++i) {
                int rwgid = b.rwgIds[i];
                sum += b.Phi[i * b.K + k] * rhsRwg[rwgid];
            }
            rhsCbf[off + k] = sum;
        }
    }
}
```

如果采用复数正交基，更严格的 Galerkin 投影通常应使用共轭转置 `Phi^H`。论文公式写成转置形式，工程实现时建议单独确认当前 EFIE/MFIE 离散的双线性形式和内积约定；如果 SVD 给出的是复正交向量，默认用 `conj(Phi)` 更稳。

### 6.5 投影近场矩阵

当前近场存储是：

```text
Z_near[rwgid] -> 若干 rwgSid 和 Zij
```

HMORM 应改成：

```text
ZnearCBF[blockA] -> 若干 near blockB 和 dense matrix Ka x Kb
```

构造方式：

```math
Z_{ab}^{CBF}
=
\Phi_a^H Z_{ab}^{RWG} \Phi_b
```

实现时不一定要显式形成完整 `Z_{ab}^{RWG}`，可以边算 RWG-RWG 近场边累加：

```cpp
for ja in 0..Ka-1:
  for kb in 0..Kb-1:
    zred = 0
    for p in rwgs(blockA):
      for q in rwgs(blockB):
        zij = computeCFIE(p, q)
        zred += conj(PhiA(p,ja)) * zij * PhiB(q,kb)
```

得到的矩阵很小，存储为 row-major：

```cpp
struct CBFNearBlock {
    int srcBlock = -1;
    int srcOffset = 0;
    int rows = 0; // K_test
    int cols = 0; // K_source
    std::vector<std::complex<double>> Z; // rows x cols
};

std::vector<std::vector<CBFNearBlock>> ZnearCBF;
```

### 6.6 预计算 CBF 远场谱

当前 `Vsmi[rwgid][kp]` 和 `Vfmj[rwgid][kp]` 已经是 RWG 级谱函数。HMORM 可以直接用它们投影得到 CBF 级谱函数：

```math
S_{b,k}(\hat{k})
=
\sum_{q\in b}
\Phi_b(q,k) Vsmi_q(\hat{k})
```

```math
V_{b,j}(\hat{k})
=
\sum_{p\in b}
\Phi_b(p,j) Vfmj_p(\hat{k})
```

代码可新增：

```cpp
void buildCBFSpectralFunctions(
    HMORMContext& ctx,
    const std::vector<std::vector<Complex3D>>& Vsmi,
    const std::vector<std::vector<Complex3D>>& Vfmj)
{
    for (auto& b : ctx.blocks) {
        int kpCount = (int)Vsmi[b.rwgIds[0]].size();
        b.S_cbf.assign(b.K, std::vector<Complex3D>(kpCount));
        b.V_cbf.assign(b.K, std::vector<Complex3D>(kpCount));

        for (int k = 0; k < b.K; ++k) {
            for (int li = 0; li < (int)b.rwgIds.size(); ++li) {
                int rwgid = b.rwgIds[li];
                auto coeff = b.Phi[li * b.K + k];
                for (int kp = 0; kp < kpCount; ++kp) {
                    b.S_cbf[k][kp].x += coeff * Vsmi[rwgid][kp].x;
                    b.S_cbf[k][kp].y += coeff * Vsmi[rwgid][kp].y;
                    b.S_cbf[k][kp].z += coeff * Vsmi[rwgid][kp].z;

                    b.V_cbf[k][kp].x += coeff * Vfmj[rwgid][kp].x;
                    b.V_cbf[k][kp].y += coeff * Vfmj[rwgid][kp].y;
                    b.V_cbf[k][kp].z += coeff * Vfmj[rwgid][kp].z;
                }
            }
        }
    }
}
```

同样，测试侧是否用共轭系数需要和内积约定一致。建议第一版保留一个开关：

```cpp
bool useConjugateProjection = true;
```

便于和原始 CFIE-MLFMA 小算例对比。

### 6.7 新增 CBF 级矩阵向量乘法

应新增一个类似 `MLFMM::matrix_solver` 的路径，例如：

```cpp
void HMORM::matrix_solver_cbf(
    int nCbf,
    std::complex<double> alpha[],
    std::complex<double> rhsCbf[],
    int itr_max,
    int mr,
    double tol_abs,
    double tol_rel);
```

它内部的 `ax` 输出 `wCbf`，流程为：

```text
alpha_CBF
  -> mexpCBF: leaf CBF source aggregation
  -> m2m: 复用当前 ProcessMLFMM::m2m
  -> m2l: 复用当前 ProcessMLFMM::m2l
  -> l2l: 复用当前 ProcessMLFMM::l2l
  -> lexpCBF: leaf CBF testing
  -> add ZnearCBF
  -> w_CBF
```

关键是新增两个函数：

```cpp
ProcessHMORM::mexpCBF(...)
ProcessHMORM::lexpCBF(...)
```

而 `m2m/m2l/l2l` 尽量复用现有 `ProcessMLFMM`。

### 6.8 CBF 级 mexp

伪代码：

```cpp
void mexpCBF(
    Sm,
    finestNodes,
    ctx,
    alpha)
{
    clear Sm leaf
    for each block b:
        int nodeIdx = b.nodeId;
        int off = ctx.cbfOffset[b.blockId];
        for k in 0..b.K-1:
            auto a = alpha[off + k];
            for kp in 0..kpCount-1:
                Sm[finest][nodeIdx][kp] += b.S_cbf[k][kp] * a;
}
```

注意这里不再访问 `node->rwgIndices` 逐 RWG 聚合。

### 6.9 CBF 级 lexp

伪代码：

```cpp
void lexpCBF(
    wCbf,
    Bm,
    ctx,
    ZnearCBF)
{
    for each block b:
        int off = ctx.cbfOffset[b.blockId];
        const auto& B = Bm[finest][b.nodeId];

        for j in 0..b.K-1:
            val = 0
            for kp in 0..kpCount-1:
                val += dot(B[kp], b.V_cbf[j][kp])
            val *= farConstant
            wCbf[off + j] += val

        for each near reduced block nb in ZnearCBF[b.blockId]:
            multiply nb.Z by alpha[srcOffset:srcOffset+cols]
            accumulate to wCbf[off:off+rows]
}
```

这里输出是 CBF 系数空间中的残差行。

### 6.10 块对角预条件器

论文的 BD 预条件器为：

```math
P_{BD}^{-1}
=
\mathrm{diag}
\left[
(Z_{11}^{CBF})^{-1},
(Z_{22}^{CBF})^{-1},
\ldots,
(Z_{MM}^{CBF})^{-1}
\right]
```

当前 `GMRES.h` 的 `mgmres` 接口只接受 `ax`，没有预条件接口。要实现右预条件，有两种方式。

第一种：新增带右预条件的 GMRES：

```text
求解 A P^{-1} y = b
最终 alpha = P^{-1} y
```

接口可为：

```cpp
template<typename AxFunc, typename ApplyRightPrecFunc>
bool mgmres_right_precond(
    int n,
    std::complex<double>* y,
    const std::complex<double>* rhs,
    int itr_max,
    int mr,
    double tol_abs,
    double tol_rel,
    AxFunc ax,
    ApplyRightPrecFunc applyPinv);
```

内部每次 Arnoldi 使用：

```cpp
z = P^{-1} v_k;
w = A z;
```

最终更新的是 `y`，求完后再：

```cpp
alpha = P^{-1} y;
```

第二种：先不改 GMRES，把左预条件作为验证版：

```text
P^{-1} A alpha = P^{-1} b
```

这实现简单，但和论文公式不是同一种右预条件。正式复现建议采用右预条件。

### 6.11 回投 RWG 电流并计算 RCS

当前 `RCSUtils::calculateRCS(...)` 需要 RWG 电流数组 `I_solve`。HMORM 求解得到的是 `alpha_CBF`，因此 RCS 前要回投：

```math
I_b^{RWG} = \Phi_b \alpha_b
```

代码可新增：

```cpp
void expandCbfSolutionToRwg(
    const HMORMContext& ctx,
    const std::complex<double>* alpha,
    std::vector<std::complex<double>>& currentRwg,
    int row)
{
    currentRwg.assign(row, 0.0);
    for (const auto& b : ctx.blocks) {
        int off = ctx.cbfOffset[b.blockId];
        for (int i = 0; i < (int)b.rwgIds.size(); ++i) {
            std::complex<double> sum = 0.0;
            for (int k = 0; k < b.K; ++k) {
                sum += b.Phi[i * b.K + k] * alpha[off + k];
            }
            currentRwg[b.rwgIds[i]] = sum;
        }
    }
}
```

然后复用当前：

```cpp
RCSUtils::calculateRCS(solver, currentRwg.data(), ...)
```

这能最大限度减少对 RCS 输出代码的扰动。

## 7. 推荐的落地路线

### 阶段 0：只做文档和开关

新增配置：

```cpp
bool useHMORM = false;
double hmormEta = 0.02;
int hmormThetaSamples = ...;
int hmormPhiSamples = ...;
bool hmormUseExtendedBlock = false;
```

先只允许：

```text
selectAlgorithm == MLFMA
integralEquation == CFIE
target == PEC
```

### 阶段 1：非扩展块 CBF 验证

目标：验证投影链路正确。

做法：

1. 每个 leaf node 作为一个 CBF block。
2. `extRwgIds = rwgIds`。
3. 生成局部 CBF。
4. 构造 `Vm_red`。
5. 构造 `ZnearCBF`。
6. 实现 `mexpCBF/lexpCBF`。
7. 求解 `alpha_CBF`。
8. 回投 RWG 电流。
9. 与当前 CFIE-MLFMA 的 RCS 和电流做对比。

这个阶段不追求论文的最大加速，只验证数学路线闭合。

### 阶段 2：周期单元分块

目标：接近论文中“最细层 cube 与 CBF block 同构”。

做法：

1. 引入阵元编号或周期单元坐标。
2. 建立 `unitCellId -> rwgIds`。
3. 让 `CBFBlock` 来自周期单元，而不是自动 octree leaf。
4. 调整或约束 octree leaf size，使一个单元尽量对应一个 leaf group。
5. 建立 `blockId -> nodeId` 映射。

如果当前输入网格没有阵元编号，可能需要从几何坐标上按周期 `dx, dy, dz` 反推单元归属。

### 阶段 3：扩展块

目标：支持电连接阵列。

做法：

1. 对每个原始块，按扩展长度 `lExt` 查找附近 RWG。
2. 构造 `extRwgIds`。
3. 解扩展块局部矩阵。
4. SVD 后用选择矩阵思想截回原始块。
5. 验证块边界电流连续性和 RCS 误差。

扩展块会引入一个实现细节：扩展块中的 RWG 可能属于多个原始块。生成 CBF 时可以共享这些 RWG，但最终 `Phi_b` 必须只覆盖原始块的 `rwgIds`，否则全局未知量会重复计数。

### 阶段 4：块对角右预条件 GMRES

目标：复现论文中 BD preconditioner 的收敛改进。

做法：

1. 对每个块保存 `ZbbCBF`。
2. 对 `ZbbCBF` 做 LU 分解或小矩阵逆。
3. 在 GMRES 中实现 `applyPinv`。
4. 求解 `A P^{-1} y = b`。
5. 输出 `alpha = P^{-1} y`。

这里不要直接对大矩阵求逆，只对每个 `K_b x K_b` 小块求逆或分解。

### 阶段 5：周期复用和哈希存储

目标：从“能跑”变成“高效”。

做法：

1. 对块间相对位移建立 key，例如 `(dx, dy, dz, typeA, typeB)`。
2. 相同 key 的近场降阶矩阵只算一次。
3. 相同几何块的 CBF 只生成一次。
4. 对 CBF 谱函数按块类型和平移相位复用。

这一步才真正接近论文中减少内存和索引开销的优势。

## 8. 与当前 PMCHWT 路径的关系

当前问题问的是 CFIE-MLFMA，因此第一版建议只做 PEC CFIE。PMCHWT 后续可以迁移，但复杂度更高。

PMCHWT 的未知量是联合向量：

```math
x =
\begin{bmatrix}
\eta J\\
M
\end{bmatrix}
```

对应代码中 `MLFMM::matrix_solver` 的 PMCHWT 分支，`x[0:N)` 和 `x[N:2N)` 分别代表两类未知量。

若要做 HMORM-PMCHWT，每个 CBF 应该作用于联合未知量：

```math
\Phi_b =
\begin{bmatrix}
\Phi_b^J\\
\Phi_b^M
\end{bmatrix}
```

近场投影也变为 2 x 2 block 投影：

```math
Z_{ab}^{CBF}
=
\Phi_a^H
\begin{bmatrix}
Z_{JJ} & Z_{JM}\\
Z_{MJ} & Z_{MM}
\end{bmatrix}
\Phi_b
```

远场也需要分别投影 `J/M` 在 `k1/k2` 下的四路谱函数。建议等 PEC CFIE-HMORM 验证完成后再做。

## 9. 最小可行实现的文件级改造建议

### 9.1 新增文件

建议新增：

```text
Code/HMORM.h
Code/HMORM.cpp
Code/CBFBlock.h
Code/ProcessHMORM.h
```

职责划分：

| 文件 | 职责 |
|---|---|
| `CBFBlock.h` | 定义 CBF 块、降阶近场块、HMORMContext |
| `HMORM.h/.cpp` | HMORM solver 主类，类似 `MLFMM` |
| `ProcessHMORM.h` | CBF 级 `mexpCBF`、`lexpCBF`、RHS 投影、电流回投 |

### 9.2 修改 `RCSExportConfig.h`

新增配置字段和校验：

```cpp
bool useHMORM = false;
double hmormEta = 0.02;
double hmormExtLengthLambda = 0.0;
int hmormMpwThetaCount = 0;
int hmormMpwPhiCount = 0;
bool hmormUseConjugateProjection = true;
```

第一版校验：

```cpp
if (useHMORM && selectAlgo != selectAlgorithm::MLFMA) return false;
if (useHMORM && integralEquation != IntegralEquation::CFIE) return false;
```

### 9.3 修改 `MLFMM.cpp`

可以不直接把 HMORM 塞进 `MLFMM::matrix_solver`，而是在 CFIE 分支中判断：

```cpp
if (cfg.useHMORM) {
    HMORM hmorm(...reuse MLFMM precomputed data...);
    hmorm.monoOrDualRcs(...);
    return;
}
```

更稳的方式是让 `HMORM` 持有或引用 `MLFMM` 中已有的：

- `octreeNodes_`
- `kp_lvl_k1`
- `WGL_k1`
- `WGL_phi_k1`
- `TF`
- `interpol_k1`
- `Vsmi`
- `Vfmj`
- `rwgs`
- `wave`
- `gausspoint`

这样第一版不重复构建 MLFMA 几何和角谱数据。

### 9.4 修改 `GMRES.h`

新增右预条件版本，而不是破坏现有 `mgmres`：

```cpp
template<typename AxFunc, typename ApplyRightPrecFunc>
bool mgmres_right_precond(...);
```

现有 CFIE-MLFMA 仍调用 `mgmres`，HMORM 调用新接口。

### 9.5 修改 RCS 计算路径

尽量复用 `RCSUtils::computeMonoStatic(...)` 的角度循环，但 HMORM 需要特殊求解流程：

```text
computeV_C -> Vm_RWG
project RHS -> Vm_CBF
solve CBF system -> alpha_CBF
expand alpha_CBF -> I_RWG
calculateRCS -> output
```

可以新增：

```cpp
RCSUtils::computeMonoStatic_HMORM(...)
RCSUtils::computeDualStatic_HMORM(...)
```

不要直接改现有 `computeMonoStatic`，以免影响已有 MLFMA、HSB、ACA-RHS 路径。

## 10. 验证方案

### 10.1 数学正确性验证

先用小模型验证：

1. 当前 CFIE-MLFMA 得到 `I_RWG_ref` 和 RCS。
2. HMORM 取非常低的 SVD 截断阈值，使 `K_b` 接近 `N_b`。
3. 回投得到 `I_RWG_hmorm`。
4. 比较：

```text
norm(I_RWG_hmorm - I_RWG_ref) / norm(I_RWG_ref)
RCS_HMORM - RCS_MLFMA
GMRES residual history
```

若 `K_b` 接近 `N_b` 时误差仍大，优先检查：

- `Phi^T` 还是 `Phi^H`；
- CBF block 的 RWG 局部顺序；
- `rwgid -> local index` 映射；
- 近场投影是否漏了自块；
- 远场常数因子是否和现有 `matrix_solver` 完全一致。

### 10.2 降阶效果验证

逐步提高 `eta0`，记录：

```text
total RWG unknowns
total CBF unknowns
compression ratio
GMRES iterations
matrix-vector time
near memory
RCS max error
RCS RMS error
```

建议先以 RCS RMS 误差小于 `0.5 dB` 或 RRMSE 小于 `5%` 作为工程阈值。

### 10.3 加速效果验证

分项计时：

```text
t_CBF_generation
t_ZnearCBF
t_SV_cbf
t_rhs_projection
t_iter
t_expand_current
t_RCS
```

注意 HMORM 有额外 CBF 生成成本。小模型可能不快，这是正常的；规模变大、多 RHS 变多时优势才会明显。

## 11. 实现优先级建议

优先级从高到低：

1. `CBFBlock` 与 `HMORMContext` 数据结构。
2. 非扩展块 CBF 生成。
3. RHS 投影与电流回投。
4. CBF 级近场矩阵投影。
5. CBF 级 `mexpCBF/lexpCBF`。
6. 无预条件 CBF GMRES 跑通。
7. BD 右预条件 GMRES。
8. 周期单元分块与块复用。
9. 扩展块。
10. PMCHWT 联合未知量扩展。

这条顺序比较稳，因为每一步都能和当前 CFIE-MLFMA 做数值对比，不会一次性把 octree、矩阵、求解器和 RCS 全部搅在一起。

## 12. 最关键的风险点

### 12.1 CBF block 与 leaf node 不一致

如果一个物理单元被拆到多个 leaf node，或者一个 leaf node 含多个物理单元，HMORM 的 CBF 复用和最细层同构映射都会变弱。快速验证版可以接受，但正式版本必须解决。

### 12.2 投影内积约定

论文公式中写 `I^T Z I`，但代码中复数离散和 SVD 正交基更常见的是 `I^H Z I`。这个差异会直接影响精度。建议实现时保留开关并通过小模型对比确定。

### 12.3 扩展块截断

扩展块生成 CBF 时包含邻域 RWG，但最终 CBF 只应落在原始块 RWG 上。若截断索引错，会出现重复自由度或边界电流不连续。

### 12.4 BD 预条件器不能直接照搬到 RWG 系统

论文 BD 预条件器是 CBF self block 的逆，不是当前 RWG 近场自块的逆。若还在 RWG 空间求解，直接套这个预条件器并不等价。

### 12.5 小规模问题可能变慢

HMORM 多了 CBF 生成、SVD、投影等预处理。只有当 `K << N`、块数大或 RHS 多时，收益才会压过预处理成本。

## 13. 总结

论文 HMORM 相较当前 CFIE-MLFMA 的本质改进是：

```text
CFIE 不变，MLFMA 层级不必推翻；
变化发生在最细层：把 RWG 自由度压缩成 CBF 自由度，
并让近场、远场、RHS、GMRES 都在 CBF 空间中工作。
```

当前代码已经有较完整的 CFIE-MLFMA 框架：

- octree 层级；
- `Vsmi/Vfmj` RWG 级谱函数；
- `TF` 转移因子；
- `mexp/m2m/m2l/l2l`；
- `Z_near` 近场；
- GMRES；
- RCS 输出。

因此实现 HMORM 时，最值得复用的是 MLFMA 的中间层传播；最需要重写的是 leaf-level 数据接口：

```text
RWG x -> CBF alpha
RWG RHS -> CBF RHS
RWG near matrix -> CBF near matrix
RWG source/testing spectrum -> CBF source/testing spectrum
CBF solution -> RWG current for RCS
```

这也是论文所谓“物理层降阶 + 算法层加速”的核心：先用物理响应子空间减少未知量，再让 MLFMA 在这个更小的空间里完成远场加速。
