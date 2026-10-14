# 基于右边向量变换的运动目标 RCS 计算方法

## 1. 研究背景

在矩量法（Method of Moments, MoM）及其快速算法 MLFMA 中，电磁散射问题离散后通常可写为线性方程组

$$
\mathbf{Z}\mathbf{I}=\mathbf{V},
$$

其中，$\mathbf{Z}$ 为阻抗矩阵，$\mathbf{I}$ 为待求表面电流系数，$\mathbf{V}$ 为入射场激励形成的右边向量。对于固定目标，$\mathbf{Z}$ 由目标几何形状、剖分网格、材料参数以及积分方程类型决定；$\mathbf{V}$ 则由入射波方向、极化方向和频率决定。

对于刚体运动目标，如果直接在每一个时刻对目标几何位置和姿态进行真实变换，则需要重复更新目标网格、重建空间层次结构、重新填充近场矩阵和远场聚合/转移相关数据。该方式计算代价很高，不适合长时间序列或高采样率的动态 RCS 计算。

本文采用的思想是：在目标做刚体运动时，不直接变换目标网格和阻抗矩阵，而是将目标运动等效为雷达入射/观测方向在目标本体系中的变化。这样可以保持阻抗矩阵 $\mathbf{Z}$ 不变，仅在每个运动时刻重新生成右边向量 $\mathbf{V}(t)$，从而实现运动目标 RCS 的快速计算。

## 2. 基本假设

该方法基于以下假设：

1. 目标为刚体运动，目标自身几何形状不随时间改变。
2. 目标材料参数和工作频率不随时间改变。
3. 每个时刻的电磁散射可看作瞬时静态散射问题。
4. 目标运动通过运动路径文件给出，路径文件中每个时刻包含目标参考点位置和姿态四元数。
5. 当前 C++ 程序中，目标运动主要通过等效改变单站入射/观测方向进入 RHS 计算；阻抗矩阵、MLFMA 层次结构和近场矩阵均保持不变。

## 3. 运动路径描述

程序使用的运动路径文件每一行格式为

```text
time_s x y z qw qx qy qz
```

其中，$t$ 为时间，$\mathbf{p}(t)=[x(t),y(t),z(t)]^T$ 为目标参考点在全局坐标系中的位置，$\mathbf{q}(t)=[q_w,q_x,q_y,q_z]$ 为目标姿态四元数。该四元数采用从目标本体系到全局坐标系的变换约定，即

$$
\mathbf{r}_g=\mathbf{R}(t)\mathbf{r}_b,
$$

其中，$\mathbf{r}_b$ 为目标本体系中的向量，$\mathbf{r}_g$ 为全局坐标系中的向量，$\mathbf{R}(t)$ 为四元数对应的旋转矩阵。

对于单位四元数

$$
\mathbf{q}=[w,x,y,z],
$$

程序采用的旋转矩阵为

$$
\mathbf{R}=
\begin{bmatrix}
1-2(y^2+z^2) & 2(xy-zw) & 2(xz+yw) \\
2(xy+zw) & 1-2(x^2+z^2) & 2(yz-xw) \\
2(xz-yw) & 2(yz+xw) & 1-2(x^2+y^2)
\end{bmatrix}.
$$

由于 $\mathbf{R}$ 为正交矩阵，因此从全局坐标系变换到目标本体系可写为

$$
\mathbf{r}_b=\mathbf{R}^T(t)\mathbf{r}_g.
$$

程序中的 `MovingTarget::rotateGlobalToBody` 正是采用 $\mathbf{R}^T$ 完成这一变换。

## 4. 单站雷达方向的等效变换

设单站雷达在全局坐标系中的位置或方向参数为

$$
\mathbf{p}_r=[x_r,y_r,z_r]^T.
$$

目标在时刻 $t_n$ 的参考点位置为

$$
\mathbf{p}_n=\mathbf{p}(t_n).
$$

则从目标指向雷达的全局视线向量为

$$
\mathbf{d}_g(t_n)=\mathbf{p}_r-\mathbf{p}_n.
$$

其距离为

$$
R_n=\|\mathbf{d}_g(t_n)\|,
$$

单位方向向量为

$$
\hat{\mathbf{s}}_g(t_n)=
\frac{\mathbf{d}_g(t_n)}{\|\mathbf{d}_g(t_n)\|}.
$$

目标姿态由旋转矩阵 $\mathbf{R}(t_n)$ 表示。为了在固定网格坐标系中计算运动目标的等效散射，需要将全局视线方向变换到目标本体系：

$$
\hat{\mathbf{s}}_b(t_n)=\mathbf{R}^T(t_n)\hat{\mathbf{s}}_g(t_n).
$$

然后根据目标本体系方向 $\hat{\mathbf{s}}_b=[s_x,s_y,s_z]^T$ 计算等效球坐标角度：

$$
\theta_b(t_n)=\arccos(s_z),
$$

$$
\phi_b(t_n)=\operatorname{atan2}(s_y,s_x).
$$

在程序中角度以度为单位输出，并将 $\phi$ 归一化到 $[0^\circ,360^\circ)$。当 $\theta=0^\circ$ 或 $\theta=180^\circ$ 时，方位角 $\phi$ 在数学上不唯一，程序将其规范化为 $0^\circ$。

因此，运动目标在时刻 $t_n$ 的单站散射计算可以等效为：固定目标在本体系中受到来自

$$
(\theta_b(t_n),\phi_b(t_n))
$$

方向的平面波照射，并在相同方向进行单站观测。

## 5. 阻抗矩阵不变的原因

MoM 离散后的阻抗矩阵元素可抽象写为

$$
Z_{mn}=
\left\langle \mathbf{f}_m,\mathcal{L}\{\mathbf{f}_n\}\right\rangle,
$$

其中，$\mathbf{f}_m$ 和 $\mathbf{f}_n$ 为 RWG 基函数，$\mathcal{L}$ 为积分算子。对于刚体目标，如果在目标本体系中观察，目标的几何拓扑、边长、三角面片面积、基函数形状、法向量相对关系均不随运动改变。因此在目标本体系中有

$$
\mathbf{Z}(t)=\mathbf{Z}_0.
$$

也就是说，只要目标不发生形变，且频率、材料和积分方程类型不变，就可以保持同一个阻抗矩阵。目标运动的影响通过等效入射方向进入右边向量：

$$
\mathbf{Z}_0 \mathbf{I}(t_n)=\mathbf{V}(t_n).
$$

这就是当前程序中“只改变右边向量，不改变阻抗矩阵”的核心依据。

## 6. 右边向量的计算公式

以 PEC 目标的 EFIE 为例，入射平面波可写为

$$
\mathbf{E}^{inc}(\mathbf{r})=
\mathbf{e}^{inc}\exp\left(-jk\hat{\mathbf{k}}\cdot\mathbf{r}\right),
$$

其中，$\hat{\mathbf{k}}$ 为入射传播方向，$\mathbf{e}^{inc}$ 为入射电场极化方向，$k$ 为波数。

对第 $m$ 个 RWG 基函数，其右边向量元素可写为

$$
V_m^{E}
=
\int_{S_m}
\mathbf{f}_m(\mathbf{r})\cdot
\mathbf{E}^{inc}(\mathbf{r})\,dS.
$$

程序中对 RWG 正负三角形分别进行高斯积分。若 $\rho_m^+(\mathbf{r})$ 和 $\rho_m^-(\mathbf{r})$ 分别表示正、负三角形上的 RWG 局部矢量函数，$l_m$ 为公共边长度，则离散计算形式可表示为

$$
V_m^{E}
\approx
l_m
\sum_{i\in T_m^+}
w_i
\left[
\mathbf{e}^{inc}\cdot\rho_m^+(\mathbf{r}_i)
\right]
\exp\left(-jk\hat{\mathbf{k}}\cdot\mathbf{r}_i\right)
+
l_m
\sum_{i\in T_m^-}
w_i
\left[
\mathbf{e}^{inc}\cdot\rho_m^-(\mathbf{r}_i)
\right]
\exp\left(-jk\hat{\mathbf{k}}\cdot\mathbf{r}_i\right).
$$

对于 MFIE 部分，程序使用入射磁场极化方向 $\mathbf{h}^{inc}$ 和面片法向量 $\mathbf{n}$ 构造磁场激励项，可抽象写为

$$
V_m^{M}
\approx
l_m
\sum_i
w_i
\left[
(\mathbf{n}\times\mathbf{h}^{inc})\cdot\rho_m(\mathbf{r}_i)
\right]
\exp\left(-jk\hat{\mathbf{k}}\cdot\mathbf{r}_i\right).
$$

CFIE 的右边向量由 EFIE 和 MFIE 加权组合：

$$
\mathbf{V}^{C}
=
\alpha\mathbf{V}^{E}
+
(1-\alpha)\mathbf{V}^{M}.
$$

因此，在运动目标计算中，不同时间 $t_n$ 的右边向量为

$$
\mathbf{V}(t_n)
=
\mathbf{V}
\left(
\theta_b(t_n),
\phi_b(t_n),
\text{pol}
\right),
$$

其中 $\text{pol}$ 表示入射极化方式。当前程序通过 `EMSource::initPW` 根据 $(\theta_b,\phi_b)$ 生成 $\hat{\mathbf{k}}$、$\mathbf{e}^{inc}$ 和 $\mathbf{h}^{inc}$，再调用 `computeV_E`、`computeV_M` 或 `computeV_C` 形成 RHS。

## 7. 准确 MLFMA-Moving 计算流程

准确计算模式下，对每一个运动时刻都重新生成 RHS，并使用同一个 MLFMA 矩阵算子求解：

$$
\mathbf{Z}_0\mathbf{I}_n=\mathbf{V}_n,
$$

其中

$$
\mathbf{V}_n=\mathbf{V}(t_n).
$$

计算流程如下：

1. 读取运动路径，得到 $\mathbf{p}_n$ 和 $\mathbf{q}_n$。
2. 计算全局视线方向 $\hat{\mathbf{s}}_g(t_n)$。
3. 通过 $\mathbf{R}^T(t_n)$ 得到本体系等效方向 $\hat{\mathbf{s}}_b(t_n)$。
4. 将 $\hat{\mathbf{s}}_b(t_n)$ 转换为 $(\theta_b,\phi_b)$。
5. 根据该方向生成右边向量 $\mathbf{V}_n$。
6. 求解 $\mathbf{Z}_0\mathbf{I}_n=\mathbf{V}_n$。
7. 用 $\mathbf{I}_n$ 计算单站 RCS。

该方法精度最高，因为每一个时刻均进行一次完整线性方程求解，但计算量也最大。

## 8. ACA-RHS-MLFMA-Moving 计算流程

对于长时间序列运动，如果时间采样点很多，逐时刻求解

$$
\mathbf{Z}_0\mathbf{I}_n=\mathbf{V}_n
$$

仍然代价较高。注意所有时刻的 RHS 可以组成矩阵

$$
\mathbf{B}=
\begin{bmatrix}
| & | & & | \\
\mathbf{V}_1 & \mathbf{V}_2 & \cdots & \mathbf{V}_{N_t} \\
| & | & & |
\end{bmatrix}
\in\mathbb{C}^{N\times N_t},
$$

其中 $N$ 为未知量个数，$N_t$ 为运动时间点数。

在许多运动目标问题中，不同时间点的 RHS 具有较强相关性，因此 $\mathbf{B}$ 可用低秩形式近似：

$$
\mathbf{B}
\approx
\mathbf{U}\mathbf{W},
$$

其中

$$
\mathbf{U}\in\mathbb{C}^{N\times r},
\qquad
\mathbf{W}\in\mathbb{C}^{r\times N_t},
\qquad
r\ll N_t.
$$

当前 ACA-RHS 算法并不先完整生成 $\mathbf{B}$，而是通过 ACA 过程按需查询某一列或某一个元素：

$$
\mathbf{B}_{:,n}=\mathbf{V}_n,
$$

$$
B_{m,n}=V_m(t_n).
$$

当 ACA 需要第 $n$ 列时，程序调用 RHS 生成函数完整计算 $\mathbf{V}_n$；当 ACA 只需要某个元素 $B_{m,n}$ 时，程序只计算对应 RWG 基函数上的 RHS 元素。最终得到低秩因子 $\mathbf{U}$ 和 $\mathbf{W}$。

随后只需要求解 $r$ 个线性方程：

$$
\mathbf{Z}_0\mathbf{X}_j=\mathbf{U}_{:,j},
\qquad j=1,2,\cdots,r.
$$

记

$$
\mathbf{X}
=
\begin{bmatrix}
| & | & & | \\
\mathbf{X}_1 & \mathbf{X}_2 & \cdots & \mathbf{X}_r \\
| & | & & |
\end{bmatrix},
$$

则第 $n$ 个时刻的电流可由低秩系数重构：

$$
\mathbf{I}_n
\approx
\sum_{j=1}^{r}
\mathbf{X}_j W_{j,n}.
$$

该方法将求解次数从 $N_t$ 次降低为 $r$ 次。当 $r\ll N_t$ 时，可以显著降低计算成本。

## 9. FSB-MLFMA-Moving 计算流程

FSB-MLFMA-Moving 的思想是利用空间带宽有限性，在角域中选取较少的准确采样方向。对运动路径中的等效角度集合

$$
\{(\theta_b(t_n),\phi_b(t_n))\}_{n=1}^{N_t}
$$

先进行角度合并和采样布局设计，只在若干代表性角度上进行准确 MLFMA 求解，然后对其他时刻的散射场进行角域插值。其基本流程为：

1. 将运动路径转换为目标本体系角度序列。
2. 合并角度上相近的运动状态。
3. 根据目标尺寸、波数和过采样因子构造 FSB 采样点。
4. 对采样点准确求解 MLFMA。
5. 对每个运动时刻插值得到对应方向的散射场。

该方法适合角域变化平滑且目标空间带宽有限的情况。

## 10. RCS 计算

在线性方程求解得到电流系数 $\mathbf{I}_n$ 后，远场散射可表示为

$$
\mathbf{E}^{sca}_\infty(t_n)
=
\mathcal{F}\{\mathbf{I}_n\},
$$

其中 $\mathcal{F}$ 表示由 RWG 电流到远场散射场的积分算子。程序中根据目标本体系的单站方向 $(\theta_b,\phi_b)$ 计算远场，并将散射场投影到垂直极化和水平极化方向，得到

$$
E_v(t_n),
\qquad
E_h(t_n).
$$

总 RCS 采用

$$
\sigma(t_n)
=
10\log_{10}
\left(
\frac{|E_v(t_n)|^2+|E_h(t_n)|^2}{4\pi}
\right),
$$

垂直极化和水平极化分量分别为

$$
\sigma_v(t_n)
=
10\log_{10}
\left(
\frac{|E_v(t_n)|^2}{4\pi}
\right),
$$

$$
\sigma_h(t_n)
=
10\log_{10}
\left(
\frac{|E_h(t_n)|^2}{4\pi}
\right).
$$

## 11. 方法特点

该方法的主要优点包括：

1. 阻抗矩阵只需构建一次，避免了每个时刻重复填充矩阵。
2. MLFMA 的层次结构、近场矩阵和转移算子可复用。
3. 目标运动通过 RHS 方向变化体现，实现较方便。
4. 可与准确逐时刻求解、ACA-RHS 低秩加速和 FSB 角域插值加速相结合。
5. 对刚体转动、平移、公转和复合姿态运动均可通过运动路径文件统一描述。

该方法的局限性包括：

1. 当前程序没有真实移动网格，因此不适用于目标发生形变的情况。
2. 当前实现主要通过等效入射/观测方向改变 RHS；如果需要严格考虑目标整体平移引入的绝对相位，需要额外确认或加入平移相位补偿项。
3. 对近场观测或距离变化显著的问题，简单方向等效可能不足，需要根据具体物理模型补充幅度和相位修正。
4. ACA-RHS 的效率依赖 RHS 矩阵的低秩性；若运动导致角域变化剧烈，低秩阶数可能升高。

## 12. 本程序中的对应关系

在当前 C++ 程序中，相关实现关系如下：

- 运动路径读取：`MovingTarget::readMovingPathFile`
- 四元数到旋转矩阵：`MovingTarget::quaternionToMatrix`
- 全局方向到本体系方向：`MovingTarget::rotateGlobalToBody`
- 方向到角度：`MovingTarget::vectorToThetaPhiDeg`
- 运动角度序列生成：`MovingTarget::computeMovingAngleSamples`
- 准确逐时刻求解：`RCSUtils::computeMonoStatic_Accurate_Moving`
- ACA-RHS 低秩加速：`RCSUtils::computeMonoStatic_ACARHS_Moving`
- FSB 角域采样加速：`RCSUtils::computeMonoStatic_FSB_Moving`

因此，当前程序实现的运动目标计算本质上是：

$$
\boxed{
\text{目标运动}
\Rightarrow
\text{等效本体系入射/观测方向}
\Rightarrow
\mathbf{V}(t)
\Rightarrow
\mathbf{Z}_0\mathbf{I}(t)=\mathbf{V}(t)
\Rightarrow
\text{RCS}(t)
}
$$

该流程避免了对目标阻抗矩阵的重复变换和重复填充，是运动目标 RCS 计算中的一种高效等效实现方式。
