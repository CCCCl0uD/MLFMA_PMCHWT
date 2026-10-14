# HMORM 公式原理拆解与应用扩展

源文献：Si Jia Guo, Yu Mao Wu, Yuan An Liu, "An Efficient Hierarchical Model Order Reduction Method for Analyzing Electromagnetic Scattering from Large-Scale Periodic Arrays", IEEE Transactions on Antennas and Propagation, DOI: 10.1109/TAP.2026.3720418.

本文提出的 HMORM（hierarchical model order reduction method）可以理解为：先用特征基函数方法（CBFM/MORM）把每个阵元或块内的大量 RWG 未知量压缩成少量宏基函数，再让 MLFMM 的最细层盒子直接对应这些 CBF 块，使迭代中的矩阵向量乘法在“降阶未知量空间”内完成。其核心收益来自两个层面：

1. 物理层降阶：用少量 CBF 表达一个单元块中的主要感应电流自由度。
2. 算法层加速：用 MLFMM 的聚合、平移、解聚合过程加速 CBF 块之间的远场耦合。

## 1. 符号与物理背景

| 符号 | 含义 |
|---|---|
| `M` | 原始块或阵元数。周期阵列中通常每个块包含一个阵元 |
| `N_m` | 第 `m` 个原始块内 RWG 基函数数目 |
| `N_{\mathrm{Ext}}^m` | 第 `m` 个扩展块内 RWG 基函数数目 |
| `K_m` | 第 `m` 个原始块保留下来的 CBF 数目，通常远小于 `N_m` |
| `N_{\mathrm{MPW}}` | 用于激励局部块的多平面波数目，`N_{\mathrm{MPW}} = 2N_\theta N_\phi` |
| `\eta_0` | SVD 截断阈值 |
| `Z_{pq}` | 第 `q` 块源电流到第 `p` 块测试函数的 RWG 阻抗子矩阵 |
| `Z_{pq}^{\mathrm{CBF}}` | CBF 降阶空间中的块间阻抗矩阵 |
| `\alpha^{\mathrm{CBF}}` | CBF 展开系数，即降阶系统中的未知量 |

论文主体针对 PEC 散射问题：闭合金属立体阵元使用 CFIE，开口金属贴片阵元使用 EFIE。作者在结论中指出，后续可扩展到复合金属-介质阵列，对应可采用 EFIE-PMCHWT 或 JMCFIE 形式。

## 2. CFIE 与 RWG 离散：原始 MoM 问题

### 2.1 CFIE 连续方程

在时间因子 `e^{-i\omega t}` 下，PEC 闭合体的 CFIE 写成

```math
-ik\eta\alpha \hat{t}\cdot\int_S \bar{G}(r,r')\cdot J(r')\,dS'
+(1-\alpha)\frac{i}{k}
\left[
\frac{1}{2}\hat{t}\cdot J(r)
-\hat{t}\cdot \hat{n}\times\nabla\times\int_S g(r,r')J(r')\,dS'
\right]
=
\alpha \hat{t}\cdot E^{inc}(r)
+(1-\alpha)\frac{i}{k}\hat{t}\cdot\hat{n}\times H^{inc}(r)
```

对应论文公式 (1)。这里第一项来自 EFIE，第二项来自 MFIE。CFIE 的意义是把 EFIE 与 MFIE 混合，避免闭合 PEC 在内谐振频率附近出现非唯一或病态问题。

格林函数为

```math
\bar{G}(r,r') =
\left(\bar{I}+\frac{\nabla\nabla}{k^2}\right)g(r,r'),
\qquad
g(r,r')=\frac{e^{ik|r-r'|}}{4\pi |r-r'|}
```

对应论文公式 (2)。

### 2.2 RWG 电流展开

MoM 中，阵列表面电流用 RWG 基函数展开：

```math
J^{\mathrm{RWG}}(r') =
\sum_{m=1}^{M}\sum_{n=1}^{N_m}
I^{\mathrm{RWG}}_{m,n} f_{m,n}(r')
```

对应论文公式 (3)。如果直接求解，未知量总数为

```math
N_{\mathrm{total}}=\sum_{m=1}^{M}N_m
```

周期阵列规模上去后，MoM 的矩阵存储和直接求解成本会迅速不可承受。HMORM 的第一步就是把每个块内的 `N_m` 个 RWG 自由度压缩为 `K_m` 个 CBF 自由度。

## 3. 扩展块 MORM：从 RWG 到 CBF

### 3.1 为什么需要扩展块

普通 CBF 如果只在原始单元块内构造，适合电分离阵元；但电连接阵列中，相邻单元共享金属边界或电流连续路径，块边界处的电流并不会自然截断。论文因此引入扩展块：对每个原始块向外扩展一圈或一定电长度，把邻域耦合对局部电流形态的影响包含进来。

作者给出的扩展长度经验范围为 `0.025\lambda` 到 `0.4\lambda`。算例中，电连接贴片阵列选择 `0.025\lambda`，相当于向外扩展一圈相邻三角面片；扩展过大反而会增加矩阵病态、引入非物理高频振荡并拉高 CPU 时间。

### 3.2 多平面波激励构造局部响应库

第 `m` 个扩展块的局部矩阵方程为

```math
Z_{\mathrm{Ext},m}^{SI} I_{\mathrm{Ext},m}^{\mathrm{MPW}}
=
V_{\mathrm{Ext},m}^{\mathrm{MPW}},
\qquad m=1,2,\ldots,N_0
```

对应论文公式 (4)。其中：

- `Z_{\mathrm{Ext},m}^{SI}` 是第 `m` 个扩展块的自阻抗矩阵；
- `V_{\mathrm{Ext},m}^{\mathrm{MPW}}` 是多平面波激励矩阵；
- `I_{\mathrm{Ext},m}^{\mathrm{MPW}}` 是这些激励下求得的局部电流响应矩阵；
- `N_0` 是需要独立构造局部块方程的块类型数。

对严格周期且几何完全相同的阵列，可以复用角块、边块、内部块等有限类型；如果块几何不同，则每个块都要独立构造，`N_0` 可能等于 `M`，降阶预处理成本明显增加。

### 3.3 SVD 提取主导电流模态

多平面波响应矩阵中包含冗余信息，论文用 SVD 压缩：

```math
I_{\mathrm{Ext},m}^{\mathrm{MPW}}
=
U_{\mathrm{Ext}}^{m} S_{\mathrm{Ext}}^{m}
\left(V_{\mathrm{Ext}}^{m}\right)^H,
\qquad m=1,2,\ldots,N_0
```

对应论文公式 (5)。保留归一化奇异值高于阈值 `\eta_0` 的左奇异向量，形成扩展块上的 CBF 电流系数矩阵 `I_{\mathrm{Ext},m}^{\mathrm{CBF}}`。

阈值选择与积分方程病态程度相关：

- 闭合 PEC 立体阵元采用 CFIE，矩阵条件较好，奇异值下降快，论文中取 `\eta_0 = 0.02`；
- 开口贴片阵元采用 EFIE，边缘奇异和病态更明显，奇异值长尾更强，论文中取 `\eta_0 = 0.001`。

### 3.4 从扩展块截回原始块

扩展块 CBF 不能直接作为最终阵列块基函数，因为扩展区域包含邻居面片。论文定义原始块和扩展块的全局 RWG 索引集合：

```math
J_\alpha^m =
\left\{\ell^{m,\alpha}_1,\ell^{m,\alpha}_2,\ldots,
\ell^{m,\alpha}_{N_\alpha^m}\right\},
\qquad \alpha\in\{\mathrm{Ext},\mathrm{Ori}\}
```

对应公式 (6)，并且有 `J_{\mathrm{Ori}}^m \subset J_{\mathrm{Ext}}^m`。

为了把全局索引映射为扩展块矩阵中的行号，定义双射：

```math
\pi^m:J_{\mathrm{Ext}}^m\rightarrow
\{1,2,\ldots,N_{\mathrm{Ext}}^m\},
\qquad
\ell^m_{\pi^m(s)}=s,\ \forall s\in J_{\mathrm{Ext}}^m
```

对应公式 (7)。

然后构造选择矩阵：

```math
C_{pq}^{m}=
\begin{cases}
1, & q=\pi^m(\ell_p^{m,\mathrm{Ori}})\\
0, & \mathrm{otherwise}
\end{cases}
```

对应公式 (8)。它的作用非常直接：从扩展块 CBF 电流中抽取属于原始块的行。

于是原始块上的 CBF 系数矩阵为

```math
I_{\mathrm{Ori},m}^{\mathrm{CBF}}
= C^m I_{\mathrm{Ext},m}^{\mathrm{CBF}},
\qquad m=1,2,\ldots,N_0
```

对应公式 (9)。这一步是扩展块技术能用于电连接阵列的关键：构造 CBF 时考虑邻域连续性，装配全局未知量时仍回到每个原始块。

## 4. 降阶矩阵方程

### 4.1 CBF 电流展开

阵列电流在 CBF 空间中展开为

```math
J^{\mathrm{CBF}}(r')=
\sum_{m=1}^{M}\sum_{k=1}^{K_m}
\alpha_{m,k}^{\mathrm{CBF}} g_{m,k}(r')
```

对应公式 (10)。其中第 `m` 块第 `k` 个 CBF 为

```math
g_{m,k}(r')=
\sum_{n=1}^{N_{\mathrm{Ori}}^m}
I_{\mathrm{Ori},m,k,n}^{\mathrm{CBF}}
f_{m,n}(r')
```

对应公式 (11)。这说明 CBF 不是新的物理网格函数，而是 RWG 基函数的局部线性组合。

### 4.2 CBF Galerkin 投影

降阶系统为

```math
Z^{\mathrm{RED}}\alpha^{\mathrm{CBF}}=V^{\mathrm{RED}}
```

对应公式 (12)。其块结构为

```math
Z^{\mathrm{RED}}=
\begin{bmatrix}
Z_{11}^{\mathrm{CBF}} & Z_{12}^{\mathrm{CBF}} & \cdots & Z_{1M}^{\mathrm{CBF}}\\
Z_{21}^{\mathrm{CBF}} & Z_{22}^{\mathrm{CBF}} & \cdots & Z_{2M}^{\mathrm{CBF}}\\
\vdots & \vdots & \ddots & \vdots\\
Z_{M1}^{\mathrm{CBF}} & Z_{M2}^{\mathrm{CBF}} & \cdots & Z_{MM}^{\mathrm{CBF}}
\end{bmatrix}
```

单个降阶块可理解为

```math
Z_{pq}^{\mathrm{CBF}}
=
\left(I_p^{\mathrm{CBF}}\right)^T
Z_{pq}
I_q^{\mathrm{CBF}},
\qquad
V_p^{\mathrm{CBF}}
=
\left(I_p^{\mathrm{CBF}}\right)^T V_p
```

对应公式 (13) 的含义。也就是说，原始 RWG 阻抗块 `Z_{pq}` 被左右投影到两个 CBF 子空间中。若 `N_p`、`N_q` 很大，而 `K_p`、`K_q` 很小，则块规模从 `N_p \times N_q` 降为 `K_p \times K_q`。

求得 `\alpha^{\mathrm{CBF}}` 后，用公式 (10)-(11) 重构电流，再由辐射积分得到散射场：

```math
E_s(r)
=
-ik\eta\int_S
\bar{G}(r,r')\cdot J^{\mathrm{CBF}}(r')\,dS'
```

对应公式 (14)。

## 5. MLFMM 中的矩阵向量乘法

### 5.1 标准 MLFMM 的近场与远场拆分

在 MLFMM 最细层，位于第 `m` 组的测试 RWG 基函数 `p` 与所有源 RWG 基函数的矩阵向量乘法可拆成近场直接项和远场多极项：

```math
\sum_{q=1}^{N_{\mathrm{total}}} Z_{pq}I_q^{\mathrm{RWG}}
=
\sum_{m'_L\in N_{m_L}}\sum_{q\in G_{m'_L}}
Z_{pq}I_q^{\mathrm{RWG}}
+
\frac{ik}{4\pi}\int
R_{m_L,p}(\hat{k})
\sum_{m'_L\notin N_{m_L}}
T_{mm'}(\hat{k}\cdot\hat{D}_{mm'})
\sum_{q\in G_{m'_L}}
F_{m'_L,q}(\hat{k})I_q^{\mathrm{RWG}}
d^2\hat{k}
```

对应公式 (15)。其中：

- `N_{m_L}` 是最细层中第 `m` 组的邻近组集合；
- 近场项直接计算并存储；
- 远场项通过聚合函数 `F`、转移函数 `T`、解聚合/接收函数 `R` 来计算。

### 5.2 RWG 聚合与解聚合函数

论文在 CFIE 框架下给出 RWG 源聚合函数与场测试函数：

```math
F_{m'_L,q}(\hat{k})
=
\int_S e^{-ik\cdot r_{qm'}}
\left(\bar{I}-\hat{k}\hat{k}\right)
\cdot f_q(r_{qm'})\,dS
```

```math
R_{m_L,p}(\hat{k})
=
\alpha\int_S e^{ik\cdot r_{pm}}
\left(\bar{I}-\hat{k}\hat{k}\right)
\cdot f_p(r_{pm})\,dS
-(1-\alpha)\hat{k}\times
\int_S e^{ik\cdot r_{pm}}
f_p(r_{pm})\times \hat{n}\,dS
```

对应公式 (16)。`F` 描述源基函数如何向外辐射平面波谱，`R` 描述测试基函数如何接收来自远场组的入射谱。

组间转移函数为

```math
T_{mm'}(\hat{k}\cdot\hat{D}_{mm'})
=
\frac{ik}{4\pi}
\sum_{l=0}^{L}
i^l(2l+1)
h_l^{(1)}(kD_{mm'})
P_l(\hat{k}\cdot\hat{D}_{mm'})
```

对应公式 (17)。这是 MLFMM 中控制两个远分离组之间传播的核心函数。

### 5.3 层级上行与下行

标准 MLFMM 的上行聚合把子层中心的 outgoing wave 平移到父层中心：

```math
F_{m'_{l-1},q}(\hat{k}_{l-1,s'})
=
e^{i\hat{k}_{l-1,s'}\cdot(r_{m'_l}-r_{m'_{l-1}})}
\sum_{s=1}^{K_l} W_{s's}F_{m'_l,q}(\hat{k}_{l,s})
```

对应公式 (18)。

远场贡献在父层上可写为

```math
I=
\int R_{m_{l-1},p}(\hat{k})\cdot B_{m_{l-1}}(\hat{k})\,d^2\hat{k}
=
\sum_{s'=1}^{K_{l-1}}
\xi_{s'}R_{m_{l-1},p}(\hat{k}_{l-1,s'})
\cdot B_{m_{l-1}}(\hat{k}_{l-1,s'})
```

对应公式 (19)。

下行解聚合把父层 incoming wave 转到子层中心：

```math
B_{m_l}(\hat{k}_{l,s})
=
\sum_{s'=1}^{K_{l-1}}
W_{s's}e^{i\hat{k}_{l-1,s'}\cdot(r_{m_l}-r_{m_{l-1}})}
B_{m_{l-1}}(\hat{k}_{l-1,s'})
\xi_{s'}/\xi_s
```

对应公式 (20)。这些是 MLFMM 的常规层级传播骨架。

## 6. HMORM 的关键嵌入：让 MLFMM 直接作用于 CBF

传统 CBFM-MLFMM 融合容易出现两套结构：一套 CBF 分块树，一套 MLFMM 八叉树。这样会引入多树构建、块-盒索引映射、矩阵装配和存储管理开销。

本文的核心设计是：把 MLFMM 最细层立方体尺寸设为阵元周期，使每个最细层 cube 对应一个 CBF block。于是：

- CBF 块分区就是 MLFMM 最细层分区；
- 每个最细层盒子包含一个单元块；
- 周期阵列中的平移不变性可直接用于复用自阻抗和近场互耦矩阵；
- MLFMM 的聚合、转移、解聚合过程不再围绕 RWG 未知量，而是围绕 CBF 展开系数。

### 6.1 CBF 源聚合函数

第 `m'` 个源组中第 `k` 个 CBF 的聚合函数为

```math
S_{m'_L,k}(\hat{k})
=
\sum_{q=1}^{N_{\mathrm{Ori}}^{m'_L}}
I_{m'_L,k,q}^{\mathrm{CBF}}
F_{m'_L,q}(\hat{k})
```

对应公式 (21)。这一步把一个 CBF 内所有 RWG 基函数的远场辐射谱线性组合起来。换句话说，MLFMM 最细层不再逐 RWG 聚合，而是逐 CBF 聚合。

### 6.2 CBF 场测试函数

第 `m` 个场组中第 `j` 个 CBF 的解聚合/接收函数为

```math
V_{m_L,j}(\hat{k})
=
\sum_{p=1}^{N_{\mathrm{Ori}}^{m_L}}
I_{m_L,j,p}^{\mathrm{CBF}}
R_{m_L,p}(\hat{k})
```

对应公式 (22)。这一步把测试侧也投影到 CBF 空间，形成 CBF 层面的接收函数。

### 6.3 CBF 空间中的近远场矩阵向量乘法

对第 `m` 个场组、第 `j` 个 CBF，降阶矩阵向量乘法为

```math
\sum_{m'_L=1}^{M}\sum_{k=1}^{K_{m'_L}}
Z_{m_L,m'_L}^{j,k,\mathrm{RED}}
\alpha_{m'_L,k}^{\mathrm{CBF}}
=
\sum_{m'_L\in N_{m_L}}\sum_{k=1}^{K_{m'_L}}
Z_{m_L,m'_L}^{j,k,\mathrm{RED}}
\alpha_{m'_L,k}^{\mathrm{CBF}}
+
\frac{ik}{4\pi}
\int
V_{m_L,j}(\hat{k})
\sum_{m'_L\notin N_{m_L}}
T_{mm'}(\hat{k}\cdot\hat{D}_{mm'})
\sum_{k=1}^{K_{m'_L}}
S_{m'_L,k}(\hat{k})
\alpha_{m'_L,k}^{\mathrm{CBF}}
d^2\hat{k}
```

对应公式 (23)。这是 HMORM 的核心公式。它和标准 MLFMM 的区别是：

| 标准 MLFMM | HMORM |
|---|---|
| 未知量是 `I_q^{RWG}` | 未知量是 `\alpha_{m,k}^{CBF}` |
| 最细层按 RWG 基函数聚合 | 最细层按 CBF 聚合 |
| 近场矩阵块为 RWG-RWG | 近场矩阵块为 CBF-CBF |
| 复杂度约 `O(MN log(MN))` | 复杂度约 `O(MK log(MK))` |

其中 `N` 是每块 RWG 数，`K` 是每块 CBF 数，通常 `K << N`。

## 7. 块对角预条件器

HMORM 使用 GMRES 求解降阶系统。为提升收敛，论文构造右预条件器：

```math
P_{\mathrm{BD}}^{-1}
=
\begin{bmatrix}
\left(Z_{11}^{\mathrm{CBF}}\right)^{-1} & 0 & \cdots & 0\\
0 & \left(Z_{22}^{\mathrm{CBF}}\right)^{-1} & \cdots & 0\\
\vdots & \vdots & \ddots & \vdots\\
0 & 0 & \cdots & \left(Z_{MM}^{\mathrm{CBF}}\right)^{-1}
\end{bmatrix}
```

对应公式 (24)。这里每个对角块是最细层自组的 CBF 自阻抗矩阵逆。它的物理含义是：先近似消除每个阵元或局部块内部的强自耦合，再由 GMRES 处理块间耦合。

预条件后的方程为

```math
Z^{\mathrm{RED}}P_{\mathrm{BD}}^{-1}
\left(P_{\mathrm{BD}}\alpha^{\mathrm{CBF}}\right)
=
V^{\mathrm{RED}}
```

对应公式 (25)。论文算例显示该预条件器对闭合 CFIE 阵列和开口 EFIE 电连接阵列都明显改善收敛。特别是 100 x 100 闭合立体阵列在 `10^{-6}` 残差下，使用 BD 预条件器为 42 次迭代，而无预条件时 503 次。

## 8. 误差指标与复杂度

RCS 误差使用 RRMSE：

```math
\delta_\Omega =
\sqrt{
\frac{
\sum_{i=1}^{N_a}
\left(\sigma_{\mathrm{cal},i}-\sigma_{\mathrm{theo},i}\right)^2
}{
\sum_{i=1}^{N_a}
\left(\sigma_{\mathrm{theo},i}\right)^2
}
},
\qquad \Omega\in\{2D,3D\}
```

对应公式 (26)。其中 `\sigma_{\mathrm{cal}}` 是 HMORM 或其他方法计算的 RCS，`\sigma_{\mathrm{theo}}` 是 FEKO 或 MLFMM 基准结果。

论文给出的复杂度对比可概括为：

```math
\mathrm{MoM}:\quad O(N^3M^3)
```

```math
\mathrm{MLFMM}:\quad O(MN\log(MN))
```

```math
\mathrm{HMORM}:\quad O(MK\log(MK))
```

更完整的 HMORM CPU 时间包含：

```math
O(9N_{\mathrm{Ext}}^3 + 9N_{\mathrm{Ext}}N_{\mathrm{MPW}}^2)
+ O(9N_{\mathrm{Ori}}^2)
+ O(9M(N_{\mathrm{Ori}}^2K+N_{\mathrm{Ori}}K^2))
+ O(MK\log(MK))
```

前几项主要是 CBF 生成、互耦和降阶矩阵预处理；当阵列规模变大时，主导项变为迭代求解中的 `O(MK log(MK))`。

## 9. 算法流程拆解

可以把 HMORM 实现成以下流水线：

1. 阵列分块：每个周期单元作为一个原始块。
2. 设定 MLFMM 最细层盒子尺寸：使最细层 cube 与 CBF block 一一对应。
3. 构造扩展块：电连接阵列向外扩展一定长度，电分离阵列可只保留原始块或更小邻域。
4. 局部多平面波激励：对典型块类型求解 `Z_{\mathrm{Ext}}^{SI}I^{MPW}=V^{MPW}`。
5. SVD 压缩：按阈值 `\eta_0` 保留主导左奇异向量。
6. 选择矩阵截断：用 `C^m` 从扩展块 CBF 中抽取原始块部分。
7. CBF 赋值复用：按角、边、内部或具体几何类型分配 CBF。
8. 近场存储：用哈希表压缩存储最细层邻近 CBF-CBF 互耦块。
9. 远场预计算：构造 CBF 聚合函数 `S` 和测试函数 `V`。
10. GMRES 求解：使用 BD 右预条件器求解降阶方程。
11. 电流重构与 RCS 计算：由 `\alpha^{CBF}` 重构 RWG 电流或直接计算散射场。

## 10. 适用边界

HMORM 最适合满足以下条件的问题：

- 几何由大量重复或近似重复单元组成；
- 每个单元的局部电流响应可由少量主模态表达；
- 块间远场耦合可由 MLFMM 有效加速；
- 需要在有限阵列中保留边缘效应，而不是假设无限周期；
- 原始 MoM 或普通 MLFMM 在未知量、内存或迭代时间上已经成为瓶颈。

需要谨慎的情况：

- 单个目标没有重复结构，CBF 复用价值下降；
- 单元尺寸过大，不能自然作为 MLFMM 最细层 cube；
- 阵元形状、姿态、材料差异很大，会使 `N_0` 接近 `M`；
- 高频强局部尖端效应或窄缝效应导致奇异值衰减慢，`K` 变大；
- 开口 EFIE 问题病态更强，需要更严格 SVD 阈值和更可靠预条件。

## 11. 除大规模有限周期阵列外的潜在应用

### 11.1 多个相同目标的 RCS 计算

这是较自然的扩展方向，但要区分两类问题。

第一类是“多个相同散射体组成的目标群”RCS。例如多个相同金属块、无人机阵列、编队目标、重复散射单元、诱饵阵列或周期/准周期排布的目标簇。若每个目标几何、姿态和网格相同，只存在平移关系，则可以把每个目标视为一个 CBF block：

- 单目标或典型目标只需生成一次 CBF；
- 目标间近场互耦可按相对位移分类复用；
- 远场耦合交给 MLFMM；
- 求解未知量从“所有目标的 RWG 总数”降为“目标数 x 每目标 CBF 数”。

这种场景和论文的有限周期阵列非常接近。若目标排列规则，收益最大；若目标位置不规则，MLFMM 仍然有效，但近场互耦矩阵的平移复用会减少。

第二类是“同一个目标在多个入射角/观察角下的 RCS 扫描”。这不是多个物理目标，而是多 RHS 问题。HMORM 仍有价值，因为 CBF、近场块、MLFMM 层级和预条件器可复用，只需改变激励向量。但论文也指出，小规模问题中直接 LU 的 MORM 在多 RHS 下可能更快；大规模问题中 HMORM 的内存优势更关键。若要高效做单站 RCS 角域扫描，建议结合 Krylov 子空间 recycling 或误差可控的角域扫描策略。

### 11.2 准周期、非周期和非正交晶格阵列

论文明确提到扩展块按块构造，因此理论上可用于 quasi-periodic 和 aperiodic arrays，例如反射阵、FSS 和非正交晶格阵列。代价是：

- 完全周期时只需少量典型块；
- 准周期或非周期时，块类型增加，CBF 生成和近场矩阵存储增加；
- 若几何变化连续，可考虑按形状参数聚类，减少独立 CBF 类型。

### 11.3 金属-介质复合阵列或介质目标群

论文结论中提出未来可用 EFIE-PMCHWT 或 JMCFIE 扩展到复合金属-介质阵列。对 PMCHWT，未知量通常同时包含等效电流 `J` 和等效磁流 `M`，因此 CBF 应作用于联合未知向量：

```math
x =
\begin{bmatrix}
J\\
M
\end{bmatrix},
\qquad
g_k^{\mathrm{CBF}} =
\sum_n I_{k,n}^{J}f_n^{J}
+
\sum_n I_{k,n}^{M}f_n^{M}
```

这对当前 PMCHWT_MLFMA 工程有直接参考价值：可以先把每个重复介质/金属-介质单元的 `[J,M]` 局部响应做 SVD，再把 MLFMA 的聚合和测试函数改写为联合电磁流的 CBF 形式。

### 11.4 大型重复结构的散射与辐射

潜在应用还包括：

- 大规模 FSS、超表面、吸波/隐身周期结构；
- 反射阵、相控阵、天线罩周期单元；
- 多端口阵列天线的有源端口响应与互耦分析；
- 地面或舰载重复结构的散射贡献快速评估；
- 由相同模块拼装的大型平台局部散射分析。

辐射问题与散射问题的不同在于 RHS 来自端口激励而不是入射平面波。论文提到可用 localized active port responses 与 mutual coupling modes 构造 radiating CBFs。

## 12. 用于“多个相同目标 RCS”的实现建议

若把 HMORM 思想迁移到多个相同目标 RCS，建议按以下方式落地：

1. 定义目标块类型：相同几何、相同姿态、相同材料、相同网格为同一类型。
2. 判断块尺度：若单个目标电尺寸较小，可把一个目标对应一个最细层 cube；若单个目标本身电大，应把目标划分为多个子块，或者采用“目标内部子块 + 目标间 MLFMM”的两级结构。
3. 构造局部 CBF：对每类目标或子块施加 MPW 激励，必要时加入近邻目标构成扩展块。
4. 建立相对位移表：规则排布时复用近场互耦；不规则排布时只复用自块和部分相同距离/姿态互耦。
5. 改写 MLFMM 最细层聚合：用公式 (21) 的 CBF 源谱替代 RWG 源谱。
6. 改写测试侧：用公式 (22) 的 CBF 接收函数替代 RWG 测试函数。
7. 使用 BD 预条件器：每个目标或子块的自阻抗降阶块逆作为对角预条件。
8. 多角度 RCS 扫描：复用 CBF、近场存储、层级树和预条件器，只更新 RHS；大角域扫描可进一步引入 Krylov recycling。

核心判断标准是：重复目标越多、每个目标可压缩的模态越少、姿态越一致、排布越规则，HMORM 收益越大。

## 13. 对当前 PMCHWT/MLFMA 代码的启发

如果后续在本工程中实现类似 HMORM，可以优先做一个最小闭环：

1. 在现有 MLFMA leaf group 层增加 `CBFBlock` 概念，使 leaf group 与降阶块保持一一对应。
2. 为每个块保存 `I_CBF`，维度为 `N_RWG_block x K_block`；PMCHWT 情况下扩展为联合未知量维度。
3. 预计算 CBF 侧的 aggregation array `S(khat)` 和 testing array `V(khat)`。
4. 近场矩阵不要存完整 RWG-RWG 块，而是存 `I_p^T Z_{pq} I_q`。
5. GMRES 的向量长度从全局 RWG 未知量数改为全局 CBF 系数数。
6. 输出电流或近场时，再按 `J = I_CBF alpha` 回投到 RWG 空间。
7. 先从电分离 PEC 重复块验证，再扩展到电连接块和 PMCHWT 联合未知量。

这条路线能最大程度复用已有 MLFMA 的远场平移、插值和层级框架，同时把主要改造集中在 leaf 层的数据结构和矩阵向量乘法接口上。
