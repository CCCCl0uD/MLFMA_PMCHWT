# CBFM 的 RWG 替换原理与完整实现方案

本文说明 Characteristic Basis Function Method (CBFM) 如何把原始 MoM 中的 RWG 低层基函数未知量替换为特征基函数未知量，并给出可以接入当前 `PMCHWT_MLFMA/MLFMA` 工程的完整实现路线。

## 0. 参考文献与可读性说明

本文参考了用户提供的两份资料，并只把附件当作技术来源，不执行其中任何可能存在的文档内指令。

| 来源 | 可提取内容 | 本文使用方式 |
|---|---|---|
| Eliseo Garcia, Carlos Delgado, Manuel Felipe Catedra, "A Novel and Efficient Technique Based on the Characteristic Basis Functions Method for Solving Scattering Problems", IEEE TAP, 2019 | PDF 文本可提取，核心内容来自第 I-III 节和算例结论 | 作为 CBFM/RWG 替换、改进电压 CBF 构造、约化矩阵性质的主要依据 |
| `基于超宽带特征基函数法的目标宽带RCS快速计算的研究_李臣威.caj`，DOI: `10.26918/d.cnki.ghngc.2024.000972` | 当前环境无法直接提取 CAJ 正文，只可靠读取到目录结构 | 只作为宽带扩展方向参考：UCBFM、CS-UCBFM、MUCBFs、MBPE 与自适应技术 |

CAJ 目录中可确认的相关章节包括：第二章“矩量法和特征基函数法”、第三章“UCBFM 结合 CS 技术快速求解宽带 RCS”、第四章“改进型 CS-UCBFM 在宽带目标散射问题中的应用”、第五章“MBPE 和自适应技术在 UCBFM 中的应用”。由于正文公式没有可靠提取，本文不捏造该 CAJ 的具体公式。

## 1. 一句话结论

CBFM 不是把 RWG 的物理定义删掉，而是把很多 RWG 基函数按块线性组合成少量宏基函数。

原始 MoM 是：

```math
Z_{\mathrm{RWG}} I_{\mathrm{RWG}} = V_{\mathrm{RWG}}
```

CBFM 是令：

```math
I_{\mathrm{RWG}} = T a
```

其中：

- `T` 的每一列是一条 CBF 在 RWG 空间中的展开系数；
- `a` 是新的未知量，即 CBF 系数；
- `T` 通常是块对角稀疏矩阵，每个块只覆盖一组局部 RWG。

代入并用高层测试函数 `W` 左投影后得到：

```math
Z_R a = V_R
```

```math
Z_R = W^H Z_{\mathrm{RWG}} T,\qquad
V_R = W^H V_{\mathrm{RWG}}
```

求出 `a` 后再回投：

```math
I_{\mathrm{RWG}} = T a
```

然后沿用现有 `RCSUtils::calculateRCS(...)` 或 `calculateRCS_PMCHWT(...)`，因为 RCS 计算仍然需要 RWG 级电流系数。

## 2. 术语表

| 符号/术语 | 含义 |
|---|---|
| `RWG` | Rao-Wilton-Glisson 低层边基函数，当前工程中对应 `RWGBase` |
| `CBF` | Characteristic Basis Function，覆盖一个块的宏基函数，是多个 RWG 的线性组合 |
| `block` | CBF 构造单元，包含一组 RWG，可对应几何分块、目标子域或 MLFMA 最细层盒子 |
| `Z_{ij}` | 第 `j` 块源 RWG 到第 `i` 块测试 RWG 的原始阻抗子矩阵 |
| `T_i` | 第 `i` 块的 CBF 展开矩阵，尺寸 `N_i x K_i` |
| `W_i` | 第 `i` 块的高层测试函数矩阵，尺寸 `N_i x K_i` |
| `K_i` | 第 `i` 块保留的 CBF 数，通常远小于块内 RWG 数 `N_i` |
| `Z_R` | CBF 约化阻抗矩阵，尺寸 `K_{\mathrm{tot}} x K_{\mathrm{tot}}` |
| `a` | CBF 系数向量，是 CBFM 真正求解的未知量 |
| `UCBFM` | Ultra-wideband CBFM，在多个频点生成或复用 CBF，用于宽带 RCS |

## 3. 原始 RWG-MoM 的形式

当前工程中，PEC EFIE/CFIE 和介质 PMCHWT 最终都离散成线性方程。以 PEC 为例，表面电流写成：

```math
J(r) = \sum_{n=1}^{N} I_n f_n(r)
```

这里 `f_n` 是第 `n` 条 RWG，`I_n` 是待求系数。Galerkin MoM 令测试函数也取 RWG，得到：

```math
\sum_{n=1}^{N} \langle L(f_n), f_m \rangle I_n
=
\langle E^{inc}, f_m \rangle
```

矩阵形式就是：

```math
Z_{\mathrm{RWG}} I_{\mathrm{RWG}} = V_{\mathrm{RWG}}
```

在当前代码中可以对应为：

| 工程位置 | 当前作用 |
|---|---|
| `Code/MyStruct.h` | `RWGBase` 保存正负三角形、公共边、自由顶点、边长等几何信息 |
| `Code/EFIE.cpp` | `EFIE::computeEFIE_Zij(...)` 计算单个 RWG-RWG EFIE 阻抗元素 |
| `Code/MFIE.cpp` | `MFIE::computeMFIE_Zij(...)` 计算单个 RWG-RWG MFIE 阻抗元素 |
| `Code/Zmartix.h` | 组装 dense MoM 矩阵或 FMM/MLFMA 近场矩阵 |
| `Code/Vm.h` | 计算 RWG 级 RHS `Vm` |
| `Code/MoM.cpp` | `matrix_solver(...)` 直接对 `Z_mom` 做矩阵向量乘 |
| `Code/MLFMM.cpp` | `matrix_solver(...)` 隐式执行 RWG 级 MLFMA 矩阵向量乘 |
| `Code/RCS.h` | 用求出的 RWG 电流计算 RCS |

CBFM 要改的关键不是 EFIE/CFIE/PMCHWT 的物理算子，而是未知量所在的线性空间。

## 4. 从 RWG 到 CBF 的数学替换

把所有 RWG 划分为 `B` 个块：

```math
\Omega = \Omega_1 \cup \Omega_2 \cup \cdots \cup \Omega_B
```

第 `i` 块包含 `N_i` 条 RWG：

```math
\{f_{i,1}, f_{i,2}, \ldots, f_{i,N_i}\}
```

第 `i` 块的第 `k` 个 CBF 定义为：

```math
g_{i,k}(r) =
\sum_{n=1}^{N_i} T_i(n,k) f_{i,n}(r),
\qquad k=1,\ldots,K_i
```

于是块内电流从：

```math
J_i(r) =
\sum_{n=1}^{N_i} I_{i,n} f_{i,n}(r)
```

改写为：

```math
J_i(r) =
\sum_{k=1}^{K_i} a_{i,k} g_{i,k}(r)
=
\sum_{n=1}^{N_i}
\left(
\sum_{k=1}^{K_i} T_i(n,k)a_{i,k}
\right) f_{i,n}(r)
```

因此：

```math
I_i = T_i a_i
```

全局拼接为：

```math
I_{\mathrm{RWG}}
=
\begin{bmatrix}
T_1 & 0 & \cdots & 0\\
0 & T_2 & \cdots & 0\\
\vdots & \vdots & \ddots & \vdots\\
0 & 0 & \cdots & T_B
\end{bmatrix}
\begin{bmatrix}
a_1\\
a_2\\
\vdots\\
a_B
\end{bmatrix}
= T a
```

如果使用同一组 CBF 做测试函数，即常规 Galerkin CBFM，则：

```math
W_i = T_i
```

约化矩阵块为：

```math
Z^R_{ij} = T_i^H Z_{ij} T_j
```

约化 RHS 为：

```math
V^R_i = T_i^H V_i
```

这就是“RWG 被替换”的本质：矩阵不再以 RWG 为行列，而是以 CBF 为行列；但是每个 CBF 内部仍然保存 RWG 展开系数。

## 5. Garcia 2019 论文中的改进 CBF 定义

Garcia 等 2019 年论文强调，可以把高层测试函数和高层展开函数分开定义，以改善约化系统条件数。

### 5.1 电压 CBF 作为测试函数

对每个孤立 block，使用一组不同方向、不同极化的平面波照射，先计算这些平面波在块内 RWG 测试函数上的 impressed voltage：

```math
V_i^{\mathrm{snap}}
=
\begin{bmatrix}
v_i^{(1)} & v_i^{(2)} & \cdots & v_i^{(P)}
\end{bmatrix}
\in \mathbb{C}^{N_i \times P}
```

其中 `P = N_{\mathrm{dir}} x N_{\mathrm{pol}}`。对 `V_i^{snap}` 做 SVD：

```math
V_i^{\mathrm{snap}} = U_i \Sigma_i Q_i^H
```

按奇异值阈值保留前 `K_i` 列：

```math
W_i = U_i(:,1:K_i)
```

这些 `W_i` 就是论文称为 voltage CBFs 的高层测试函数。它们是在低层测试函数系数空间中正交归一的：

```math
W_i^H W_i = I
```

### 5.2 由电压 CBF 生成电流 CBF

然后对同一个孤立 block 的自阻抗矩阵 `Z_{ii}`，逐列求解：

```math
Z_{ii} J_{i,k} = W_{i,k}
```

矩阵形式：

```math
Z_{ii} T_i = W_i
```

因此：

```math
T_i = Z_{ii}^{-1} W_i
```

`T_i` 的每一列就是一个电流 CBF：它表示当块受到第 `k` 个 voltage CBF 激励时，块内低层 RWG 上诱导出来的电流分布。

### 5.3 为什么自块约化矩阵变成单位阵

对同一块 `i=j`，改进方法的约化自块为：

```math
Z^R_{ii} = W_i^H Z_{ii} T_i
```

代入：

```math
Z_{ii}T_i = W_i
```

得到：

```math
Z^R_{ii} = W_i^H W_i = I
```

这就是 Garcia 2019 论文的关键优势：最强的块内自耦合项被精确正交化，约化矩阵的对角块天然是单位阵。论文算例显示，这能减少显著互耦项数量、降低迭代次数，并减少自块耦合项预处理时间。

### 5.4 常规 CBFM 与改进 CBFM 的区别

| 项目 | 常规 CBFM | Garcia 2019 改进 CBFM |
|---|---|---|
| 测试函数 | 通常用电流 CBF，`W_i=T_i` | 用 SVD 后的 voltage CBF，`W_i=U_i` |
| 展开函数 | 由局部激励响应/SVD 得到 | 解 `Z_{ii}T_i=W_i` 得到电流 CBF |
| 自块矩阵 | `T_i^H Z_{ii} T_i`，不是单位阵 | `W_i^H Z_{ii}T_i=I` |
| 预处理 | 需要计算自块约化矩阵 | 自块可直接填单位阵 |
| 收敛 | 通常优于 RWG-MoM | 对病态或大尺寸问题更稳 |

## 6. CBFM 完整算法流程

### 6.1 输入

实现 CBFM 至少需要：

- 全局 `rwgs`；
- 原始积分方程类型：EFIE、CFIE 或 PMCHWT；
- 原始算子参数：`k1`、`eta1`、介质情况下的 `k2`、`eta2`；
- 高斯积分点；
- 分块策略；
- 局部平面波激励采样；
- SVD 截断阈值；
- 求解器参数。

### 6.2 分块

分块目标是让每个 block 的电尺寸不太小也不太大。

建议第一版采用以下策略：

```text
若接入 dense MoM:
    可按几何包围盒或 OCTree leaf 合并成 block

若接入 MLFMA:
    优先复用 octreeNodes_[maxLevel_] 作为 block
    或者把若干相邻 leaf 合成一个 CBF block
```

每个 block 保存：

```cpp
struct CBFBlock {
    int blockId = -1;
    std::vector<int> rwgIds;
    std::vector<int> reducedIds;
    Point center{};
    double radius = 0.0;
};
```

`reducedIds` 是该 block 内 CBF 在全局约化向量 `a` 中的编号。

### 6.3 生成局部激励快照

对第 `i` 个 block，构造多个入射方向和极化：

```text
for theta in thetaSamples:
    for phi in phiSamples:
        for pol in {theta-pol, phi-pol}:
            compute local RHS v_i(theta, phi, pol)
```

形成：

```math
V_i^{snap} =
\left[
v_i^{(1)}, v_i^{(2)}, \ldots, v_i^{(P)}
\right]
```

对 PEC EFIE 可以直接复用 `RHS::computeV_E_entry(...)`。对 CFIE 复用 `RHS::computeV_C_entry(...)`。对 PMCHWT 则生成 `2N_i` 长度的联合 RHS：

```math
v_i =
\begin{bmatrix}
v_i^J\\
v_i^M
\end{bmatrix}
```

### 6.4 SVD 截断

对每个 block 的快照矩阵做薄 SVD：

```math
V_i^{snap}=U_i\Sigma_iQ_i^H
```

选择：

```math
K_i = \max \left\{
k \mid \sigma_k / \sigma_1 \ge \tau_{\mathrm{svd}}
\right\}
```

建议参数：

| 场景 | `tau_svd` 初值 |
|---|---|
| PEC CFIE | `1e-2` 到 `2e-2` |
| PEC EFIE | `1e-3` 到 `1e-2` |
| PMCHWT | `1e-3` 起步 |
| 宽带 UCBFM | `1e-3`，再用 RCS 误差回调 |

### 6.5 生成 CBF 展开矩阵

有两种可选实现。

方案 A：常规 CBFM。

```math
T_i = U_i(:,1:K_i),\qquad W_i=T_i
```

方案 B：Garcia 2019 改进 CBFM，建议作为目标实现。

```math
W_i = U_i(:,1:K_i)
```

```math
Z_{ii}T_i=W_i
```

如果 `Z_{ii}` 是 `N_i x N_i`，需要对同一个 `Z_{ii}` 做一次 LU/QR 分解，然后解 `K_i` 个 RHS。

### 6.6 构造约化 RHS

原始 RHS 仍先按 RWG 计算：

```math
V_{\mathrm{RWG}}
```

然后分块投影：

```math
V_i^R = W_i^H V_i
```

全局拼接：

```math
V_R =
\begin{bmatrix}
V_1^R\\
V_2^R\\
\vdots\\
V_B^R
\end{bmatrix}
```

### 6.7 构造约化阻抗矩阵

对任意 block 对 `(i,j)`：

```math
Z^R_{ij}=W_i^H Z_{ij} T_j
```

如果 `i==j` 且使用 Garcia 改进 CBFM：

```math
Z^R_{ii}=I
```

无需再计算。

如果 `i!=j`，可采用直接矩阵块投影：

```text
Zij_raw = compute raw RWG-RWG coupling between block i and block j
tmp = Zij_raw * T_j
Zij_red = W_i^H * tmp
```

### 6.8 求解和回投

解约化系统：

```math
Z_R a = V_R
```

然后回投 RWG 电流：

```math
I_i = T_i a_i
```

拼成全局：

```math
I_{\mathrm{RWG}} = T a
```

对 PMCHWT：

```math
\begin{bmatrix}
I_J\\
I_M
\end{bmatrix}
=
T a
```

再分别传给：

```cpp
RCSUtils::calculateRCS_PMCHWT(
    solver,
    I_J,
    I_M,
    RCS_pre,
    solver.wave.k1(),
    solver.wave.k2(),
    solver.wave.eta1(),
    solver.wave.eta2(),
    theta,
    phi);
```

## 7. 当前工程的 dense MoM 第一版实现

建议先做 dense MoM 版，因为它最容易验证“RWG 到 CBF”的替换是否正确。

### 7.1 新增文件建议

```text
Code/CBFM.h
Code/CBFM.cpp
```

核心数据结构：

```cpp
namespace CBFM {

struct Options {
    double svdTol = 1.0e-3;
    int thetaSamples = 6;
    int phiSamples = 12;
    int minCbfPerBlock = 1;
    int maxCbfPerBlock = 64;
    bool useVoltageTesting = true;
};

struct Basis {
    std::vector<CBFBlock> blocks;
    std::vector<std::vector<std::complex<double>>> T; // N_global x K_total, sparse-first later
    std::vector<std::vector<std::complex<double>>> W; // N_global x K_total
    int originalUnknowns = 0;
    int reducedUnknowns = 0;
};

} // namespace CBFM
```

第一版为了清楚可以先用 dense `T/W`，验证通过后再改成 block-local 存储。

### 7.2 建议接口

```cpp
namespace CBFM {

Basis buildBasisPEC(
    const std::vector<RWGBase>& rwgs,
    const std::vector<CBFBlock>& blocks,
    const EMSource& wave,
    const gaussPoints& gp,
    int integralEquation,
    const Options& opt);

void buildReducedSystemDense(
    const std::vector<std::vector<std::complex<double>>>& Z,
    const std::vector<std::complex<double>>& rhs,
    const Basis& basis,
    std::vector<std::vector<std::complex<double>>>& Zred,
    std::vector<std::complex<double>>& rhsRed);

void expandCurrent(
    const Basis& basis,
    const std::vector<std::complex<double>>& alpha,
    std::vector<std::complex<double>>& currentRwg);

} // namespace CBFM
```

### 7.3 在 `MoM` 中接入

当前 `MoM` 构造流程是：

```text
compute Z_mom
compute Vm for each incident angle
matrix_solver(row, I_solve, Vm, ...)
calculateRCS(I_solve)
```

CBFM 版改成：

```text
compute Z_mom
build CBF basis T/W
for each incident angle:
    compute Vm at RWG level
    project rhsRed = W^H Vm
    solve Zred alpha = rhsRed
    expand I_solve = T alpha
    calculateRCS(I_solve)
```

其中 `Zred` 如果频率和几何不变，只需构造一次。单站 RCS 角扫只变化 RHS，可以复用 `Zred`。

### 7.4 dense 版伪代码

```cpp
CBFM::Basis basis = CBFM::buildBasisPEC(
    rwgs, blocks, wave, gausspoint, integralEquType_, opt);

std::vector<std::vector<std::complex<double>>> Zred;
std::vector<std::complex<double>> rhsRed;

CBFM::buildReducedSystemDense(Z_mom, Vm, basis, Zred, rhsRed);

std::vector<std::complex<double>> alpha(basis.reducedUnknowns, 0.0);
solveDenseOrGMRES(Zred, rhsRed, alpha);

std::vector<std::complex<double>> I_solve;
CBFM::expandCurrent(basis, alpha, I_solve);

RCSUtils::calculateRCS(
    *this,
    I_solve.data(),
    RCS_pre,
    wave.k1(),
    wave.eta1(),
    theta,
    phi);
```

### 7.5 dense 版验证标准

第一阶段用小目标验证：

1. 直接 MoM 求 `I_{\mathrm{RWG}}^{MoM}`。
2. CBFM 求 `I_{\mathrm{RWG}}^{CBFM}=T a`。
3. 比较：

```math
\epsilon_I =
\frac{\|I_{\mathrm{RWG}}^{MoM}-I_{\mathrm{RWG}}^{CBFM}\|_2}
{\|I_{\mathrm{RWG}}^{MoM}\|_2}
```

4. 比较 RCS：

```math
\epsilon_{\mathrm{RCS}} =
\sqrt{
\frac{\sum_p |\sigma_p^{MoM}-\sigma_p^{CBFM}|^2}
{\sum_p |\sigma_p^{MoM}|^2}
}
```

5. 检查维度压缩比：

```math
\rho = K_{\mathrm{tot}} / N_{\mathrm{RWG}}
```

如果 `rho` 太大，说明分块过小、SVD 阈值过低，或激励快照冗余不足。

## 8. MLFMA-CBFM 实现路线

dense CBFM 的瓶颈是仍然需要 `Z_mom`。大规模问题应把约化矩阵向量乘写成隐式 MLFMA。

### 8.1 最直接但较慢的隐式实现

最小改造可以这样做：

```math
y = Z_R a = W^H Z_{\mathrm{MLFMA}} T a
```

每次 GMRES 的 `ax_red(a)`：

```text
1. expand x_rwg = T a
2. call existing MLFMA RWG matvec: y_rwg = Z x_rwg
3. project y_red = W^H y_rwg
```

优点：几乎不用改现有 `ProcessMLFMM::mexp/m2m/m2l/l2l`。

缺点：每次迭代仍在 RWG 级做 MLFMA 聚合和测试，未知量减少了，但 leaf 层工作量没有完全减少。

这个版本适合作为 MLFMA-CBFM 的正确性桥接。

### 8.2 真正高效的 leaf-level CBFM

高效版本要把 MLFMA 最细层的源聚合和测试函数都改成 CBF 级。

当前 RWG 级 `mexp` 是：

```cpp
Sm[node][kp] += Vsmi[rwgid][kp] * x[rwgid];
```

CBF 级应预计算：

```math
S_{i,k}(\hat{k}) =
\sum_{n=1}^{N_i} T_i(n,k) F_{i,n}(\hat{k})
```

代码概念：

```cpp
Sm[node][kp] += CBFSourceSpectrum[cbfId][kp] * alpha[cbfId];
```

当前 RWG 级测试是：

```cpp
val += dot(B[node][kp], Vfmj[rwgid][kp]);
```

CBF 级应预计算：

```math
R_{i,k}(\hat{k}) =
\sum_{n=1}^{N_i} \overline{W_i(n,k)} R_{i,n}(\hat{k})
```

代码概念：

```cpp
val += dot(B[node][kp], CBFTestingSpectrum[cbfId][kp]);
```

这样 `matrix_solver` 的向量长度就是 `K_total`，leaf 聚合和 leaf 测试也都是 CBF 数，而不是 RWG 数。

### 8.3 CBF 近场矩阵

当前近场矩阵是 RWG CSR：

```text
Z_near[rwgid][j]
Z_near_id[rwgid][j]
```

CBF 级近场需要变成：

```math
Z^{R,near}_{ij}=W_i^H Z^{near}_{ij} T_j
```

存储建议：

```cpp
struct CBFNearBlock {
    int testBlock = -1;
    int sourceBlock = -1;
    int rowOffset = 0;
    int colOffset = 0;
    int rowCount = 0;
    int colCount = 0;
    std::vector<std::complex<double>> values; // row-major K_i x K_j
};
```

矩阵向量乘：

```cpp
for (const CBFNearBlock& nb : nearBlocks) {
    for (int r = 0; r < nb.rowCount; ++r) {
        std::complex<double> sum = 0.0;
        for (int c = 0; c < nb.colCount; ++c) {
            sum += nb.values[r * nb.colCount + c] *
                   alpha[nb.colOffset + c];
        }
        y[nb.rowOffset + r] += sum;
    }
}
```

如果使用 Garcia 改进 CBFM，`testBlock == sourceBlock` 的自块可直接填单位阵：

```cpp
if (i == j && opt.useVoltageTesting) {
    ZijRed = identity(K_i);
}
```

### 8.4 MLFMA-CBFM 的 `ax_red`

最终目标形式：

```cpp
auto ax_red = [this](int kTotal,
                    const std::complex<double>* alpha,
                    std::complex<double>* y) {
    fill(y, y + kTotal, std::complex<double>(0.0, 0.0));

    ProcessMLFMM::clearSmBm(Sm, Bm);
    ProcessCBFMLFMM::mexpCBF(Sm, cbfBlocks, cbfSourceSpectrum, alpha);
    ProcessMLFMM::m2m(Sm, octreeNodes_, kp_lvl_k1, interpol_k1,
                      maxLevel_, levelSpan, wave.k1());
    ProcessMLFMM::m2l(Sm, Bm, octreeNodes_, TF, kp_lvl_k1,
                      WGL_k1, WGL_phi_k1, phi_level_k1,
                      maxLevel_, levelSpan);
    ProcessMLFMM::l2l(Bm, octreeNodes_, kp_lvl_k1, phi_level_k1,
                      interpol_k1, maxLevel_, levelSpan, wave.k1());
    ProcessCBFMLFMM::lexpCBF(y, cbfBlocks, cbfTestingSpectrum, Bm);
    ProcessCBFMLFMM::addNearCBF(y, nearBlocks, alpha);
};
```

## 9. PMCHWT 中的 CBFM

当前工程支持介质 PMCHWT，未知量是联合向量：

```math
x =
\begin{bmatrix}
J\\
M
\end{bmatrix}
```

其中当前代码使用：

```text
x[0:N)   = electric current coefficients
x[N:2N) = magnetic current coefficients
```

PMCHWT 的 CBF 不应只压缩 `J` 或只压缩 `M`，而应压缩联合未知量。

第 `i` 块的联合 CBF：

```math
g_{i,k}^{PMCHWT} =
\begin{bmatrix}
g_{i,k}^{J}\\
g_{i,k}^{M}
\end{bmatrix}
=
\sum_{n=1}^{N_i}
\begin{bmatrix}
T_i^J(n,k)\\
T_i^M(n,k)
\end{bmatrix}
f_{i,n}
```

局部快照和自阻抗为：

```math
V_i^{snap} \in \mathbb{C}^{2N_i \times P}
```

```math
Z_{ii}^{PMCHWT} \in \mathbb{C}^{2N_i \times 2N_i}
```

Garcia 改进式生成：

```math
W_i = U_i(:,1:K_i)
```

```math
Z_{ii}^{PMCHWT} T_i = W_i
```

约化块：

```math
Z_{ij}^R = W_i^H Z_{ij}^{PMCHWT} T_j
```

回投后拆分：

```cpp
std::vector<std::complex<double>> xRwg2N;
CBFM::expandCurrent(basis, alpha, xRwg2N);

const std::complex<double>* I_J = xRwg2N.data();
const std::complex<double>* I_M = xRwg2N.data() + solver.row;
```

## 10. UCBFM 宽带扩展

CAJ 题名和目录表明其重点是“超宽带特征基函数法的目标宽带 RCS 快速计算”，并进一步结合 CS、MUCBFs、MBPE 和自适应技术。由于正文未可靠提取，下面给出符合 UCBFM 通用思路、且适合当前工程的实现路线。

### 10.1 多频点 CBF 库

给定频带：

```math
f \in [f_{\min}, f_{\max}]
```

选择少量训练频点：

```math
\{f_1,f_2,\ldots,f_L\}
```

对每个 block 和每个训练频点生成局部电流 CBF 快照：

```math
T_i^{snap} =
\left[
T_i(f_1), T_i(f_2), \ldots, T_i(f_L)
\right]
```

再做一次 SVD：

```math
T_i^{snap}=U_i\Sigma_iQ_i^H
```

保留主导列作为 UCBF：

```math
T_i^{UCBF}=U_i(:,1:K_i^{U})
```

### 10.2 每个频点求解

对目标频点 `f`：

```text
1. 更新 k/eta/介质参数
2. 构造或隐式计算 Z(f)
3. 投影 Z_R(f)=W^H Z(f)T
4. 投影 V_R(f)=W^H V(f)
5. 求解 alpha(f)
6. 回投 I_RWG(f)=T alpha(f)
7. 计算 RCS(f)
```

### 10.3 与 CS/MBPE 的关系

可以把 UCBFM 作为物理降阶层，再把 CS 或 MBPE 作为频率/角度采样层：

- UCBFM 减少每个频点每个角度的未知量；
- CS 用少量频点或角度样本恢复稀疏表示下的宽带响应；
- MBPE 用有理函数或 Padé 型模型从少量频点外推局部频率响应；
- 自适应策略用于决定何时补充训练频点或观测角点。

第一版不建议同时实现 UCBFM、CS 和 MBPE。推荐顺序是：

```text
single-frequency CBFM
-> multi-RHS CBFM angular sweep
-> UCBFM multi-frequency basis
-> adaptive frequency sampling
-> CS/MBPE acceleration
```

## 11. 与现有代码的改造清单

### 11.1 最小可验证版本

| 文件 | 修改 |
|---|---|
| `Code/CBFM.h` | 新增 CBFM 数据结构和接口 |
| `Code/CBFM.cpp` | 实现分块、快照生成、SVD、约化、回投 |
| `Code/RCSExportConfig.h` | 增加 `useCBFM`、`cbfSvdTol`、`cbfThetaSamples`、`cbfPhiSamples`、`cbfBlockMode` |
| `Code/MoM.cpp` | 在 dense MoM 路径中加入 CBFM 分支 |
| `Code/MLFMM.cpp` | 后续加入 CBFM-MLFMA 分支 |
| `Code/MLFMA.vcxproj` | 加入新源文件 |
| `Code/MLFMA.vcxproj.filters` | 加入工程筛选器 |

### 11.2 需要的线性代数能力

当前工程没有明显的 SVD/LU 库封装。建议二选一：

1. 引入 Eigen：
   - `Eigen::BDCSVD` 或 `JacobiSVD`；
   - `Eigen::PartialPivLU` / `ColPivHouseholderQR`；
   - 开发最快。

2. 使用 LAPACK/MKL：
   - `zgesdd` / `zgesvd`；
   - `zgetrf` + `zgetrs`；
   - 性能更强，但工程配置更重。

第一版建议 Eigen，因为 CBFM 的主要风险在算法结构，不在线性代数库本身。

### 11.3 配置建议

```cpp
bool useCBFM = false;
bool cbfUseVoltageTesting = true;
double cbfSvdTol = 1.0e-3;
int cbfThetaSamples = 6;
int cbfPhiSamples = 12;
int cbfMinPerBlock = 1;
int cbfMaxPerBlock = 64;
std::string cbfBlockMode = "octree_leaf";
bool exportCBFInfo = true;
bool exportCBFCurrent = false;
```

## 12. 关键实现细节

### 12.1 复数内积必须使用共轭

约化投影应使用 Hermitian transpose：

```math
W^H Z T
```

也就是代码中：

```cpp
sum += std::conj(W[row][k]) * value;
```

如果误用普通转置，非实数问题会产生系统性相位错误。

### 12.2 CBF 列归一化

建议保存 `T_i` 前对列做归一化，或至少记录列范数。若使用 Garcia 改进式，`W_i` 已正交归一，`T_i=Z_{ii}^{-1}W_i` 不一定单位范数。求解稳定性和 RCS 回投都必须使用同一套 `T_i`，不要在约化矩阵和回投阶段使用不同归一化。

### 12.3 block 边界与扩展块

对电流在 block 边界不连续敏感的问题，普通 block CBF 可能误差偏大。可借鉴 HMORM/扩展块思想：

```text
生成 CBF 时使用 extended block
装配全局 T 时只截回 original block 的 RWG 行
```

这样能让局部 CBF 感知邻域耦合，同时避免全局未知量重复覆盖。

### 12.4 CBF 数量下限

每块至少保留一个 CBF：

```cpp
Ki = std::clamp(Ki, opt.minCbfPerBlock, opt.maxCbfPerBlock);
```

但如果某块 `N_i` 很少，必须保证：

```math
K_i \le N_i
```

### 12.5 RHS 快照方向采样

平面波方向要覆盖单位球。第一版可以使用简单网格：

```text
theta = pi * (t + 0.5) / Ntheta
phi   = 2*pi * p / Nphi
pol   = theta-pol and phi-pol
```

更好的做法是使用 Fibonacci sphere 或 Lebedev 点，使方向分布更均匀。

### 12.6 多 RHS 复用

单站 RCS 扫描中，几何、频率、阻抗矩阵不变，只有入射方向变化。CBFM 的收益很明显：

```text
Z_R 只构造一次
每个角度只计算 V_R、解 alpha、回投电流、算 RCS
```

如果使用直接分解 `Z_R`，多角度扫描会很快。

## 13. 推荐实现阶段

### 阶段 1：dense PEC EFIE/CFIE CBFM

目标：

- 支持 `MoM + PEC + EFIE/CFIE`；
- 使用 dense `Z_mom`；
- 生成 `T/W`；
- 构造 `Z_R`；
- 求解、回投、RCS；
- 和原始 MoM 对比。

验收：

```text
小球/小 almond:
    CBFM RCS 与 MoM 接近
    K_total < N_RWG
    Garcia 改进式 Z_R 自块接近单位阵
```

### 阶段 2：dense PMCHWT CBFM

目标：

- 支持 `MoM + PMCHWT`；
- 联合压缩 `[J; M]`；
- 回投后调用 `calculateRCS_PMCHWT`。

验收：

```text
介质小球:
    与原始 PMCHWT MoM 对比 RCS
    检查 J/M 两部分没有错位
```

### 阶段 3：MLFMA-CBFM 桥接版

目标：

- GMRES 在 CBF 空间求解；
- 每次 matvec 内部执行 `x_rwg=T*a`、现有 MLFMA matvec、`y=W^H*y_rwg`。

验收：

```text
结果与 dense CBFM 一致
不要求 leaf 层加速最大化
```

### 阶段 4：MLFMA-CBFM 高效版

目标：

- 预计算 CBF 源谱和测试谱；
- 近场存储 CBF-CBF 块；
- leaf 层循环从 RWG 改成 CBF；
- GMRES 全程只见 `K_total`。

验收：

```text
内存低于 RWG-MLFMA
迭代 matvec 时间低于桥接版
RCS 与原始 MLFMA 或 dense CBFM 接近
```

### 阶段 5：UCBFM 宽带版

目标：

- 多训练频点生成 UCBF；
- 对频带内目标频点复用 UCBF；
- 输出宽带 RCS；
- 后续再接入自适应采样、CS 或 MBPE。

## 14. 最容易出错的地方

| 风险 | 表现 | 排查 |
|---|---|---|
| 忘记共轭转置 | RCS 相位/幅度异常 | 检查所有 `W^H` 是否用了 `conj` |
| block 索引错位 | 回投电流乱序 | 导出 `rwgId -> block/local/globalReducedId` |
| PMCHWT J/M 错位 | 介质 RCS 完全不对 | 明确 `0:N` 和 `N:2N` 的映射 |
| CBF 阈值太松 | 误差大 | 降低 `svdTol`，增加快照方向 |
| CBF 阈值太紧 | 降阶收益小 | 提高 `svdTol` 或增大 block |
| 分块过小 | CBF 数接近 RWG 数 | 合并 leaf 或按电尺寸分块 |
| 分块过大 | 局部矩阵求解昂贵 | 限制每块 RWG 上限 |
| 宽带 UCBF 训练频点不足 | 频带边缘误差大 | 自适应补点 |

## 15. 一个清晰的最小闭环

如果马上开始写代码，建议第一轮只做下面这个闭环：

```text
输入: 当前 MoM PEC CFIE 单站 RCS 小算例

1. 正常构造 Z_mom
2. 用 OCTree leaf 或简单几何划分生成 blocks
3. 每个 block 采样 6 x 12 x 2 个平面波 RHS
4. SVD 得到 W_i
5. 解 Z_ii T_i = W_i 得到 T_i
6. 拼装 Z_R = W^H Z_mom T
7. 每个扫描角:
       compute Vm
       rhs_R = W^H Vm
       solve Z_R alpha = rhs_R
       I = T alpha
       calculateRCS(I)
8. 输出:
       CBF 数量
       每块 Ki
       Z_R 自块误差 ||Z_ii^R-I||
       RCS 对比误差
```

这个闭环一旦跑通，就已经证明“RWG 替换为 CBF”的核心机制正确。随后再把 dense `Z_mom` 替换成 MLFMA 隐式 `Z*x`，就是工程规模化的问题。

## 16. 和现有 HMORM 文档的关系

项目中已有 `docs/HMORM公式原理拆解与应用扩展.md`，其中 HMORM 可以看成 CBFM 在周期/阵列场景下与 MLFMM 深度融合的版本。

本文更基础，关注：

- CBFM 本身如何从 RWG 系数空间降阶；
- Garcia 2019 的 voltage-CBF 改进为什么能让自块变成单位阵；
- 如何先在当前 MoM/MLFMA 代码中实现一个可验证闭环；
- 如何扩展到 PMCHWT 和 UCBFM 宽带 RCS。

因此建议实际开发顺序是先实现本文的 dense CBFM，再复用 HMORM 文档中的 leaf-level CBF-MLFMA 思路做高效化。
