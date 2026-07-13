# MLFMM_PEC 与 Fortran MLFMA 对比文档

## 背景

本次对比目标是排查当前 C++ `MLFMM` 使用 CFIE 计算 PEC 目标结果不正确的问题。你提到近期参考了：

```text
D:\MyCode\Integral_Solver_V3.2_OMP
```

中的 Fortran 项目，对 MLFMM 的多极子截断数和 `d2k` 权重进行了修改。因此本文重点比较：

1. 多极子截断数 `L` 与角向采样数。
2. `d2k = WGL[t] * WGL_phi` 权重所在位置。
3. PEC-CFIE 的聚合/转移/配置/解聚合流程。
4. CFIE 的 EFIE/MFIE 组合方式。

主要参考文件：

- Fortran: `D:\MyCode\Integral_Solver_V3.2_OMP\Integral_Solver_V3.2_OMP\Integral Solver Kernel\MLFMA\MLFMA.f90`
- C++: `D:\MyCode\PMCHWT_MLFMA\MLFMA\Code\AlgoMLFMM.h`
- C++: `D:\MyCode\PMCHWT_MLFMA\MLFMA\Code\ProcessMLFMM.h`
- C++: `D:\MyCode\PMCHWT_MLFMA\MLFMA\Code\MLFMM.cpp`
- C++: `D:\MyCode\PMCHWT_MLFMA\MLFMA\Code\ProcessFMM.h`

## 结论摘要

当前 C++ 的 MLFMM_PEC 与 Fortran 参考实现至少存在两个高优先级差异，足以导致 CFIE 结果错误。

### 高优先级问题 1：角向采样数与 Fortran 不一致

Fortran 中：

```fortran
num_multipole_model = ceiling(k*dl + 6.0 * (k*dl)**(1.0/3.0))
num_theta_multipole = num_multipole_model + 1
num_phai_multipole  = 2 * num_multipole_model + 1
```

当前 C++ 中：

```cpp
const int L = computeMultipoleOrder(cubeLen, wave.k1_abs());
const int thetaCount = L;
const int phiCount = 2 * L;
```

这比 Fortran 少了：

- `theta`: 少 1 个点
- `phi`: 少 1 个点

如果 `L` 表示 Fortran 的 `num_multipole_model`，那么 C++ 应使用：

```cpp
const int thetaCount = L + 1;
const int phiCount = 2 * L + 1;
```

当前 `L, 2L` 会改变：

- 最细层 `Vsmi/Vfmj` 的维度。
- 每一层 `kp_lvl` 的维度。
- `WGL` 和 `WGL_phi`。
- 插值矩阵。
- 转移矩阵。
- 最终角谱积分。

这是最可能导致 PEC-CFIE 远场错误的原因之一。

### 高优先级问题 2：`d2k` 权重位置与 Fortran 不等价

Fortran 的 `ML` 函数中，M2L 完成后立即对本层局部展开乘积分权重：

```fortran
B(1,:,:,i_cube)=B(1,:,:,i_cube)+transf(:,:,index_transfer)*C(1,:,:,index_cube)
B(2,:,:,i_cube)=B(2,:,:,i_cube)+transf(:,:,index_transfer)*C(2,:,:,index_cube)

B(1,:,:,i_cube)=B(1,:,:,i_cube)*weight_integral(:,:)
B(2,:,:,i_cube)=B(2,:,:,i_cube)*weight_integral(:,:)
```

当前 C++ 的 `m2l` 不乘权重：

```cpp
Bmvec[kp_i].x += TFvec[kp_i] * Smvec[kp_i].x;
Bmvec[kp_i].y += TFvec[kp_i] * Smvec[kp_i].y;
Bmvec[kp_i].z += TFvec[kp_i] * Smvec[kp_i].z;
```

而是在 `lexpPEC` 最后统一乘最细层权重：

```cpp
const double d2k = WGL[tt] * WGL_phi;
val += d2k * lexpKpDot(B[tt * pN + pp], V[tt * pN + pp]);
```

这与 Fortran 的逐层 `ML` 乘权重不等价。

原因是 MLFMM 中每一层都有自己的角向网格和权重。粗层的 M2L 贡献在 `L2L` 中会被反插值到更细层。如果粗层贡献没有在粗层先乘自己的 `weight_integral`，而是等到最细层再乘最细层 `d2k`，则粗层贡献使用了错误的积分权重。

因此，若目标是严格对齐 Fortran，C++ 应恢复为：

```cpp
Bm[level][node][kp] += TF[level][vec][kp] * Sm[level][dr][kp] * weight(level, kp);
```

其中：

```cpp
int t = kp_i / pN_level;
double weight = WGL[rev][t] * WGL_phi[rev];
```

并且 `lexpPEC` 中不再重复乘 `d2k`。

## 详细对比

## 1. 多极子截断数 `L`

### Fortran

Fortran 在 `MLFMA.f90` 中定义：

```fortran
real(kind=DP),parameter::dig_multipole=6.0E0_DP
```

截断数计算：

```fortran
dl = sqrt(3.0E0_DP) * dd
k = abs(k0)
num_multipole_model = ceiling(k*dl + dig_multipole * (k*dl)**(1.0E0_DP/3.0E0_DP))
```

### C++

当前 C++：

```cpp
inline int computeMultipoleOrder(double cubeLen, double kAbs) {
    const double kD = kAbs * std::sqrt(3.0) * cubeLen;
    return static_cast<int>(std::ceil(kD + 6.0 * std::cbrt(kD)));
}
```

### 对比结论

`L` 的公式本身一致：

```text
L = ceil(kD + 6 * cbrt(kD))
```

这里不是主要问题。

真正的问题在于下一步角向采样数。

## 2. `theta/phi` 采样数

### Fortran

PEC-CFIE 分支：

```fortran
this%octree(i)%num_theta_multipole = this%octree(i)%num_multipole_model + 1
this%octree(i)%num_phai_multipole  = 2 * this%octree(i)%num_multipole_model + 1
```

### C++

当前 C++：

```cpp
const int thetaCount = L;
const int phiCount = 2 * L;
```

### 差异

| 项目 | Fortran | 当前 C++ | 是否一致 |
|---|---:|---:|---|
| 截断数 | `L` | `L` | 一致 |
| theta 点数 | `L + 1` | `L` | 不一致 |
| phi 点数 | `2L + 1` | `2L` | 不一致 |
| phi 步长 | `2pi / (2L + 1)` | `2pi / (2L)` | 不一致 |

### 建议

如果继续以 Fortran 为基准，应改回：

```cpp
const int thetaCount = L + 1;
const int phiCount = 2 * L + 1;
```

并同步检查 PMCHWT 的 k2 分支也使用同样规则。

## 3. 球面采样方向和权重

### Fortran

Fortran 的 `cal_sph_k`：

```fortran
call fn0(fn, ak, num_theta_multipole)
dphai = pi2 / num_phai_multipole
phai = (i_phai - 0.5E0_DP) * dphai

sph_k(1,i_siti,i_phai) = sqrt(1 - fn(i_siti)**2) * cos(phai)
sph_k(2,i_siti,i_phai) = sqrt(1 - fn(i_siti)**2) * sin(phai)
sph_k(3,i_siti,i_phai) = -fn(i_siti)

weight_integral(i_siti,i_phai) = ak(i_siti) * dphai
```

### C++

当前 C++：

```cpp
my_math::gauss_legendre(thetaCount, XGL_temp, WGL_temp);
WGL_k1.push_back(WGL_temp);

for (int j = thetaCount - 1; j >= 0; --j) {
    theta.push_back(std::acos(XGL_temp[j]));
}

const double dPhi = 2.0 * Pi / static_cast<double>(phiCount);
WGL_phi_k1.push_back(dPhi);
phi.push_back((j + 0.5) * dPhi);

kp = {
    sin(theta[t]) * cos(phi[p]),
    sin(theta[t]) * sin(phi[p]),
    cos(theta[t])
};
```

### 对比结论

采样形式大体一致，都是：

```text
phi = (j + 0.5) * dPhi
d2k = GaussWeight(theta) * dPhi
```

但因为 C++ 的 `thetaCount/phiCount` 与 Fortran 不一致，所以实际采样点和权重仍然不一致。

## 4. `d2k` 权重位置

### Fortran

Fortran 在每一层 `ML` 里乘权重：

```fortran
level%B = ML(..., level%transf, level%C, level%weight_integral)
```

`ML` 内部：

```fortran
B = B + transf * C
B = B * weight_integral
```

之后 `LL` 使用已经带权重的 `B` 做向下配置：

```fortran
level%A = LL(..., level_1%B or level_1%A, ..., B_child, k1)
```

最终 `Lexp` 只是求和：

```fortran
res(edgeID) = sum(REC(:,:,:,edgeID) * A(:,:,:,i_cube))
```

### C++

当前 C++：

```cpp
m2l(...);  // no weight
l2l(...);
lexpPEC(... WGL_k1[0], WGL_phi_k1[0] ...);  // finest-level weight only
```

### 对比结论

不一致，并且不是简单的代码位置差异，而是数学意义不同。

Fortran 是：

```text
每一层 M2L 贡献使用该层自己的 quadrature weight
```

当前 C++ 是：

```text
所有层的 M2L 贡献经过 L2L 后，统一使用最细层 quadrature weight
```

这会导致粗层贡献权重错误。

### 建议修正

若要对齐 Fortran，建议：

1. `m2l` 恢复接收：

```cpp
const std::vector<std::vector<double>>& WGL,
const std::vector<double>& WGL_phi,
const std::vector<std::vector<double>>& phi_level
```

2. `m2l` 内部恢复逐层权重：

```cpp
const int pN = static_cast<int>(phi_level[rev].size());
for (int kp_i = 0; kp_i < kps; ++kp_i) {
    const int t = kp_i / pN;
    const double d2k = WGL[rev][t] * WGL_phi[rev];
    Bmvec[kp_i].x += TFvec[kp_i] * Smvec[kp_i].x * d2k;
    Bmvec[kp_i].y += TFvec[kp_i] * Smvec[kp_i].y * d2k;
    Bmvec[kp_i].z += TFvec[kp_i] * Smvec[kp_i].z * d2k;
}
```

3. `lexpPEC` 改回不乘 `d2k`：

```cpp
val += lexpKpDot(B[kpIdx], V[kpIdx]);
```

4. `MLFMM.cpp` 的 PEC 调用改回：

```cpp
ProcessMLFMM::m2l(Sm, Bm, octreeNodes_, TF,
    kp_lvl_k1, WGL_k1, WGL_phi_k1, phi_level_k1, maxlvl, levelSpan);
```

5. `lexpPEC` 调用不再传 `WGL/WGL_phi`。

## 5. Mexp 聚合

### Fortran

```fortran
C(:,:,:,i_cube)=C(:,:,:,i_cube)+RAD(:,:,:,EdgeID)*vecx(EdgeID)
```

### C++

```cpp
Sm[maxlvl][nodeNum][ki].x += Vsmi[id][ki].x * x[xOffset + id];
Sm[maxlvl][nodeNum][ki].y += Vsmi[id][ki].y * x[xOffset + id];
Sm[maxlvl][nodeNum][ki].z += Vsmi[id][ki].z * x[xOffset + id];
```

### 对比结论

流程一致：最细层按 RWG 汇聚源展开。

但数据表示不同：

- Fortran 使用 `theta/phi` 两个切向分量。
- C++ 使用 `x/y/z` 三个笛卡尔分量。

这可以等价，但要求 `Vsmi/Vfmj` 中的投影关系与 Fortran 的 `sph_theta/sph_phi` 完全一致。

## 6. MM 上行聚合

### Fortran

```fortran
rho(:)=cube_child(i_cube)%c(:)-cube_father(parentID)%c(:)
expkr(:,:)=exp(-cj*k*kr(:,:))
C_father = C_father + interpolation(C_child) * expkr
```

### C++

```cpp
double vec2F[3] = {
    node->center.x - son->center.x,
    node->center.y - son->center.y,
    node->center.z - son->center.z
};
std::complex<double> e_ = exp(-J * k_wave * dot(vec2F, kp));
```

### 对比结论

这里需要谨慎核对符号。

Fortran 使用：

```text
rho = child - father
exp(-j k khat · rho)
```

C++ 使用：

```text
vec2F = father - child
exp(-j k khat · vec2F)
```

即：

```text
exp(+j k khat · (child - father))
```

这看起来与 Fortran 相反。

不过 C++ 的 `Vsmi` 初始相位也与 Fortran 的 `RAD` 表达可能存在相反定义，因此不能只凭这一处断定错误。但如果当前结果不正确，建议把 MM/L2L 的相位符号列为第二优先级核查项。

## 7. ML/M2L 转移

### Fortran

Fortran `cal_transf_kernel` 对同一个转移向量同时生成两个方向：

```fortran
value(i,j,1) = tmp_value
value(i,j,2) = tmp_value1
```

其中 `tmp_value1` 使用 `f(-k_rho_dot)`，对应正反方向复用。

`ML` 使用：

```fortran
index_transfer = cube(i_cube)%index_transf(i_far_cube)
B = B + transf(:,:,index_transfer) * C(:,:,index_cube)
```

### C++

C++ 的 `OCTree::getDRVec()` 直接把带符号的：

```cpp
node->center - dr->center
```

作为独立向量存入 `octreeNodesDRvec_`，`TF[level][vecid]` 只存一个方向。

### 对比结论

两种实现可以等价：

- Fortran: 一个几何向量存正反两个 transfer。
- C++: 正反方向作为两个 signed vector 分别存。

但必须确认 `dRVecId` 的向量方向与 `TF` 中使用的 `rho` 方向和 `m2l` 中的 `drid` 一致。

建议增加一次 debug 输出，随机选一个 node/dr：

```text
node.center - dr.center
octreeNodesDRvec_[level][dRVecId]
```

确认二者完全一致。

## 8. LL 下行配置

### Fortran

```fortran
rho(:)=cube_child(i_cube)%c(:)-cube_father(parentID)%c(:)
expkr(:,:)=exp(cj*k*kr(:,:))
A = B_child + anterpolation(B_father_or_A_father * expkr)
```

### C++

```cpp
r = child.center - father.center
e_ = exp(-J * k_wave * dot(kp, r));
dst += e_ * Bm2Vec * interpolationWeight;
```

### 对比结论

同样存在符号疑点。

Fortran 是：

```text
exp(+j k khat · (child - father))
```

C++ 是：

```text
exp(-j k khat · (child - father))
```

这和 MM 的情况相反。若 C++ 使用的 RAD/REC 相位约定不同，可能整体抵消；但如果要严格对齐 Fortran，需要系统检查：

1. `Vsmi` 的初始相位。
2. `Vfmj` 的测试端相位。
3. `m2m` 符号。
4. `l2l` 符号。
5. transfer kernel 使用 hankel 类型和 `j^l` 因子。

不要只单独改一个符号。

## 9. Lexp 解聚合

### Fortran

```fortran
res(edgeID)=sum(REC(:,:,:,edgeID)*A(:,:,:,i_cube))
res(:)=res(:)*etar0/(4.0E0_DP*lambda*lambda)
```

因为：

```text
lambda = 2*pi/k
etar0 / (4 lambda^2) = eta * k^2 / (16 pi^2)
```

### C++

```cpp
val += d2k * lexpKpDot(B[kpIdx], V[kpIdx]);
val *= (k_ * eta_ / (4.0 * Pi)) * (k_ / (4.0 * Pi));
```

系数部分：

```text
(k eta / 4pi) * (k / 4pi) = eta * k^2 / (16 pi^2)
```

### 对比结论

外部系数一致。

问题仍然是：

```text
d2k 不应只在最细层 Lexp 乘
```

如果恢复 Fortran 的 M2L 层级权重，`lexpPEC` 中应移除 `d2k`。

## 10. PEC-CFIE 的 RAD/REC 组合

### Fortran

Fortran 中：

```fortran
this%RAD = RAD
this%REC = conjg(RAD)

if (.not. justEFIE) then
    this%REC = alpha_CFIE * this%REC + (1.0 - alpha_CFIE) * REC_MFIE
end if
```

即：

```text
source side: RAD is EFIE aggregation
test side: REC = alpha * EFIE_REC + (1-alpha) * MFIE_REC
```

### C++

当前 C++：

```cpp
Vsmi = KernelVsmi_Efie
Vfmj = KernelVfmj_Cfie
```

`KernelVfmj_Cfie` 中：

```cpp
res = alpha * efie + (1 - alpha) * mfie;
```

### 对比结论

结构上与 Fortran 一致：

- 源端仍是 EFIE 聚合。
- 测试端混合 CFIE。

但需要继续核对 MFIE 项的符号。Fortran 中 `REC_MFIE` 的核心是：

```fortran
tmp_REC = tmp_REC - cross(conjg(Ctemp), tri%n)
tmp_REC = cross(k_MLFMA, tmp_REC) * lm
```

C++ 中：

```cpp
cross(n_cross_pf, n_, p_f);
cross(k_cross_n_cross_pf, kp, n_cross_pf);
mfie = k_cross_n_cross_pf;
```

这里是否等价取决于：

- `p_f` 正负面元的符号处理。
- `normalP/normalN` 方向。
- 是否需要 `conj`。
- `khat` 使用方向。

如果修复采样数和层级权重后 CFIE 仍不对，再重点查这一项。

## 推荐修复顺序

### 第一步：恢复 Fortran 的角向采样数

在 `Code/AlgoMLFMM.h` 中把 k1/k2 两处：

```cpp
const int thetaCount = L;
const int phiCount = 2 * L;
```

改为：

```cpp
const int thetaCount = L + 1;
const int phiCount = 2 * L + 1;
```

这是最直接、最确定的不一致。

### 第二步：恢复 Fortran 的 M2L 层级权重

在 `Code/ProcessMLFMM.h` 中让 `m2l` 重新按层接收并乘：

```cpp
WGL[rev][t] * WGL_phi[rev]
```

并从 `lexpPEC` 中移除最细层 `d2k`。

这一点非常重要，因为粗层 M2L 贡献必须使用粗层权重。

### 第三步：重新编译并只测 PEC-EFIE

先把 `alpha = 1` 或选择 EFIE 路径跑通。EFIE 只涉及 EFIE REC，可避免 MFIE 符号干扰。

如果 EFIE 恢复正确，再测 CFIE。

### 第四步：若 CFIE 仍错误，再查 MFIE REC 符号

重点比较：

- Fortran `tmp_REC = tmp_REC - cross(conjg(Ctemp), tri%n)`
- Fortran `tmp_REC = cross(k_MLFMA, tmp_REC) * lm`
- C++ `KernelVfmj_Cfie` 的 `cross(n, p_f)` 和 `cross(k, n_cross_pf)` 顺序

### 第五步：若 EFIE 仍错误，再系统查相位符号

不要单独改 `m2m` 或 `l2l`。应同时列出：

| 阶段 | Fortran | C++ |
|---|---|---|
| RAD/Vsmi | `exp(-j k khat · (r-center))` | 当前等价表达需核对 |
| MM/m2m | `exp(-j k khat · (child-father))` | 当前看似相反 |
| LL/l2l | `exp(+j k khat · (child-father))` | 当前看似相反 |
| Lexp/Vfmj | `conjg(RAD)` 或 CFIE REC | 当前 kernel 表达需核对 |

只有整体相位约定确定后再改符号。

## 建议的最小代码修复方向

如果目标是先尽快恢复与 Fortran 一致，建议本轮只做两个改动。

### 1. `AlgoMLFMM.h`

```cpp
const int thetaCount = L + 1;
const int phiCount = 2 * L + 1;
```

k1 和 k2 两处都改。

### 2. `ProcessMLFMM.h`

恢复 `m2l` 权重：

```cpp
inline void m2l(...,
    const std::vector<std::vector<kp_Point>>& kp_lvl,
    const std::vector<std::vector<double>>& WGL,
    const std::vector<double>& WGL_phi,
    const std::vector<std::vector<double>>& phi_level,
    int maxLevel_, int levelSpan)
{
    ...
    const int pN = static_cast<int>(phi_level[rev].size());
    ...
    const int t = kp_i / pN;
    const double d2k = WGL[rev][t] * WGL_phi[rev];
    Bmvec[kp_i].x += TFvec[kp_i] * Smvec[kp_i].x * d2k;
    Bmvec[kp_i].y += TFvec[kp_i] * Smvec[kp_i].y * d2k;
    Bmvec[kp_i].z += TFvec[kp_i] * Smvec[kp_i].z * d2k;
}
```

同时 `lexpPEC` 改为：

```cpp
val += lexpKpDot(B[kpIdx], V[kpIdx]);
```

## 最终判断

当前 C++ MLFMM_PEC 与 Fortran 的主要偏差不是截断公式本身，而是：

1. `L` 之后的角向点数少了 `+1`。
2. `d2k` 从 Fortran 的“每层 M2L 后乘”改成了“最细层 Lexp 统一乘”。

这两点会直接影响 PEC-CFIE 的远场矩阵向量乘结果。建议优先恢复这两处，再继续检查 CFIE 的 MFIE 测试端符号。
