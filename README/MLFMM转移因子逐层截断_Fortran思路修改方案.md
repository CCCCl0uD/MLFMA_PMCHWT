# MLFMM 转移因子逐层截断修改方案

## 1. 问题背景

当前 C++ 代码中，`AlgoFMM::computeBase()` 先用最细层盒子边长计算一个全局的 `L_k1` / `L_k2`：

```cpp
const double highLevelLen = octreeNodes[maxLevel_][0]->cubeLength;
double D_ = std::sqrt(3.0) * highLevelLen;
L_k1 = ceil(k1 * D_ + 6 * cbrt(k1 * D_));
L_k2 = ceil(k2 * D_ + 6 * cbrt(k2 * D_));
```

随后各层角向采样点数通过倍增得到：

```cpp
int n = (1 << i) * L_k1;
phiCount = 2 * n;
```

但是在 `AlgoMLFMM::computeTransferFactorMatrix()` 中，转移因子展开阶数仍然使用同一个 `L`：

```cpp
ComplexBessel::spherical_hankel_2_array(L, kR, h2_array);
for (int l = 0; l <= L; ++l) {
    ...
}
```

这意味着：虽然 C++ 的粗层角向采样点数变多了，但转移核的球汉克尔/勒让德展开仍然只用最细层的截断阶数。对低频、PEC 或浅层问题，这个误差可能不明显；对高频介质 PMCHWT，特别是 `k2 = k1 * sqrt(eps_r)` 较大时，粗层转移核会被低阶截断，容易导致 MLFMM 的矩阵向量乘法误差增大，最终表现为 GMRES/CGS 不收敛或残差停滞。

## 2. Fortran 代码思路

Fortran 中 `find_multipole_num_PMCHW()` 对每一层分别计算多极子模式数：

```fortran
do i=1,this%max_level
    ! 1空间
    this%octree(i)%num_multipole_model =
        cal_multipole_num(this%octree(i)%cube_edge_dl,k1)

    this%octree(i)%num_theta_multipole =
        this%octree(i)%num_multipole_model+1

    this%octree(i)%num_phai_multipole =
        2*this%octree(i)%num_multipole_model+1

    ! 2空间
    this%octree(i)%num_multipole_model2 =
        cal_multipole_num(this%octree(i)%cube_edge_dl,k2)

    this%octree(i)%num_theta_multipole2 =
        this%octree(i)%num_multipole_model2+1

    this%octree(i)%num_phai_multipole2 =
        2*this%octree(i)%num_multipole_model2+1
end do
```

模式数公式为：

```fortran
dl = sqrt(3.0E0_DP) * dd
num_multipole_model =
    ceiling(k*dl + dig_multipole*(k*dl)**(1.0E0_DP/3.0E0_DP))
```

其中：

- `dd` 是当前层盒子边长。
- `dl` 是当前层盒子对角线长度。
- `k` 是当前空间波数的模，PMCHWT 中有 `k1` 和 `k2` 两套。
- `dig_multipole` 在当前 C++ 中等价于硬编码的 `6.0`。

转移矩阵生成时，Fortran 在 `set_transf_PMCHW()` 中按当前层取尺寸：

```fortran
do i_level=3,this%max_level
    num_theta_multipole = this%octree(i_level)%num_theta_multipole
    num_phai_multipole  = this%octree(i_level)%num_phai_multipole
    allocate(my_transf(num_theta_multipole,num_phai_multipole,transf_size))

    num_theta_multipole2 = this%octree(i_level)%num_theta_multipole2
    num_phai_multipole2  = this%octree(i_level)%num_phai_multipole2
    allocate(my_transf2(num_theta_multipole2,num_phai_multipole2,transf_size))

    my_transf(:,:,num_transf:num_transf+1) =
        cal_transf_kernel(this,k1,num_theta_multipole,num_phai_multipole,
                          this%octree(i_level)%sph_k,rho(:))

    my_transf2(:,:,num_transf:num_transf+1) =
        cal_transf_kernel(this,k2,num_theta_multipole2,num_phai_multipole2,
                          this%octree(i_level)%sph_k2,rho(:))
end do
```

因此 Fortran 的链路是：

```text
当前层盒长
  -> 当前层 L
  -> 当前层 theta/phi 点数
  -> 当前层 kp 网格
  -> 当前层转移矩阵尺寸
  -> 当前层转移因子截断阶数
```

## 3. C++ 修改目标

C++ 应改成保存每一层的多极子模式数，而不是只保存一个 `L_k1/L_k2`。

建议新增两个数组：

```cpp
std::vector<int> L_level_k1;
std::vector<int> L_level_k2;
```

层序建议与当前 `kp_lvl_k1/kp_lvl_k2` 保持一致：

```text
index 0              -> maxLevel_ 最细层
index 1              -> maxLevel_ - 1
...
index levelSpan - 1  -> 第 2 层
```

这样现有 `rev = levelSpan - 1 - level` 的索引方式可以继续使用。

## 4. 修改 `MLFMM.h` 和 `FMM.h`

在 `MLFMM` 类中保留 `L_k1/L_k2` 作为兼容字段，但新增每层模式数：

```cpp
int L_k1, L_k2;
std::vector<int> L_level_k1;
std::vector<int> L_level_k2;
```

如果 `FMM` 也共用 `AlgoFMM::computeBase()`，则 `FMM.h` 中也要同步增加：

```cpp
std::vector<int> L_level_k1;
std::vector<int> L_level_k2;
```

建议含义：

- `L_k1/L_k2`：最细层模式数，仅用于兼容现有输出或 FMM 单层路径。
- `L_level_k1/L_level_k2`：MLFMM 每一层真实模式数。

## 5. 修改 `AlgoFMM::computeBase()` 接口

当前接口只返回 `L_k1/L_k2`。建议增加两个输出参数：

```cpp
std::vector<int>& L_level_k1,
std::vector<int>& L_level_k2
```

修改后的函数签名建议为：

```cpp
inline void computeBase(
    const std::vector<std::vector<OCTree::Node*>>& octreeNodes,
    const EMSource& wave,
    const int integralEquType_,
    int maxLevel_,
    int rwgSize,
    int& L_k1,
    int& L_k2,
    int& row,
    int& levelSpan,
    std::vector<int>& L_level_k1,
    std::vector<int>& L_level_k2,
    std::vector<std::vector<double>>& WGL_k1,
    std::vector<std::vector<double>>& WGL_k2,
    std::vector<double>& WGL_phi_k1,
    std::vector<double>& WGL_phi_k2,
    std::vector<std::vector<double>>& theta_level_k1,
    std::vector<std::vector<double>>& theta_level_k2,
    std::vector<std::vector<double>>& phi_level_k1,
    std::vector<std::vector<double>>& phi_level_k2,
    std::vector<std::vector<kp_Point>>& kp_lvl_k1,
    std::vector<std::vector<kp_Point>>& kp_lvl_k2)
```

## 6. 在 `computeBase()` 中逐层计算 L

新增一个小工具函数，避免 `k1/k2` 重复写：

```cpp
inline int computeMultipoleOrder(double cubeLen, double kAbs) {
    const double kd = kAbs * std::sqrt(3.0) * cubeLen;
    return static_cast<int>(std::ceil(kd + 6.0 * std::cbrt(kd)));
}
```

然后按当前 C++ 的层序，从最细层到第 2 层填充：

```cpp
levelSpan = maxLevel_ - 1;
L_level_k1.clear();
L_level_k2.clear();

for (int i = 0; i < levelSpan; ++i) {
    const int treeLevel = maxLevel_ - i;
    const double cubeLen = octreeNodes[treeLevel][0]->cubeLength;
    const int L = computeMultipoleOrder(cubeLen, wave.k1_abs());
    L_level_k1.push_back(L);
}

L_k1 = L_level_k1.empty() ? 0 : L_level_k1[0];
```

PMCHWT 时再填 `k2`：

```cpp
if (integralEquType_ == 2) {
    for (int i = 0; i < levelSpan; ++i) {
        const int treeLevel = maxLevel_ - i;
        const double cubeLen = octreeNodes[treeLevel][0]->cubeLength;
        const int L = computeMultipoleOrder(cubeLen, wave.k2_abs());
        L_level_k2.push_back(L);
    }
    L_k2 = L_level_k2.empty() ? 0 : L_level_k2[0];
}
```

## 7. 角向采样点数改成 Fortran 规则

Fortran 是：

```text
thetaCount = L + 1
phiCount   = 2 * L + 1
```

C++ 建议也改成这个规则。对 `k1`：

```cpp
for (int i = 0; i < levelSpan; ++i) {
    const int L = L_level_k1[i];
    const int thetaCount = L + 1;
    const int phiCount = 2 * L + 1;

    AlgoFMM::gauss_legendre(thetaCount, XGL_temp, WGL_temp);
    WGL_k1.push_back(WGL_temp);

    theta.clear();
    for (int j = thetaCount - 1; j >= 0; --j) {
        theta.push_back(std::acos(XGL_temp[j]));
    }
    theta_level_k1.push_back(theta);

    phi.clear();
    const double dPhi = 2.0 * Pi / static_cast<double>(phiCount);
    WGL_phi_k1.push_back(dPhi);
    for (int j = 0; j < phiCount; ++j) {
        phi.push_back((j + 0.5) * dPhi);
    }
    phi_level_k1.push_back(phi);

    std::vector<kp_Point> current_lvl_kp;
    current_lvl_kp.reserve(thetaCount * phiCount);
    for (int t = 0; t < thetaCount; ++t) {
        for (int p = 0; p < phiCount; ++p) {
            current_lvl_kp.emplace_back(
                std::sin(theta[t]) * std::cos(phi[p]),
                std::sin(theta[t]) * std::sin(phi[p]),
                std::cos(theta[t]));
        }
    }
    kp_lvl_k1.push_back(current_lvl_kp);
}
```

`k2` 完全同理，只是使用 `L_level_k2`、`WGL_k2`、`theta_level_k2`、`phi_level_k2`、`kp_lvl_k2`。

注意：这一步会改变 `kp_lvl` 每层大小，插值矩阵和工作区大小会自动随 `kp_lvl` 改变；但所有依赖 `theta_level/phi_level` 的地方都要确保没有假设 `phiCount == 2 * thetaCount`。

## 8. 修改转移因子函数接口

当前接口：

```cpp
AlgoMLFMM::computeTransferFactorMatrix(
    TF, octreeNodesDRvec_, kp_lvl_k1, maxLevel_, L_k1, wave.k1());
```

建议改成：

```cpp
AlgoMLFMM::computeTransferFactorMatrix(
    TF, octreeNodesDRvec_, kp_lvl_k1, L_level_k1, maxLevel_, wave.k1());
```

函数签名改为：

```cpp
inline void computeTransferFactorMatrix(
    std::vector<std::vector<std::vector<std::complex<double>>>>& TF,
    const std::vector<std::vector<OCTree::nodePoint>>& octreeNodesDRvec,
    const std::vector<std::vector<kp_Point>>& kp_lvl,
    const std::vector<int>& L_level,
    int maxLevel_,
    std::complex<double> k_wavenumber)
```

在函数内部，每层取自己的 `L`：

```cpp
for (int level = 0; level < maxLevel_ - 1; ++level) {
    const int rev = span - 1 - level;
    const int L = L_level[rev];
    ...
    ComplexBessel::spherical_hankel_2_array(L, kR, h2_array);

    for (int l = 0; l <= L; ++l) {
        ...
    }
}
```

这样第 2 层用第 2 层的粗盒子模式数，第 3 层用第 3 层模式数，最细层用最细层模式数，与 Fortran 思路一致。

## 9. 修改 `MLFMM.cpp` 调用点

构造函数中 `computeBase()` 调用要加入两个数组：

```cpp
AlgoFMM::computeBase(
    octreeNodes_, wave, integralEquType_, maxLevel_,
    static_cast<int>(rwgs.size()),
    L_k1, L_k2, row, levelSpan,
    L_level_k1, L_level_k2,
    WGL_k1, WGL_k2, WGL_phi_k1, WGL_phi_k2,
    theta_level_k1, theta_level_k2,
    phi_level_k1, phi_level_k2,
    kp_lvl_k1, kp_lvl_k2);
```

PEC/CFIE/EFIE 的 `k1` 转移矩阵：

```cpp
AlgoMLFMM::computeTransferFactorMatrix(
    TF, octreeNodesDRvec_, kp_lvl_k1, L_level_k1,
    maxLevel_, wave.k1());
```

PMCHWT 的 `k1/k2` 转移矩阵：

```cpp
AlgoMLFMM::computeTransferFactorMatrix(
    TF, octreeNodesDRvec_, kp_lvl_k1, L_level_k1,
    maxLevel_, wave.k1());

AlgoMLFMM::computeTransferFactorMatrix(
    TF2, octreeNodesDRvec_, kp_lvl_k2, L_level_k2,
    maxLevel_, wave.k2());
```

## 10. 修改 `FMM.cpp` 调用点

如果 `FMM` 仍然共用 `AlgoFMM::computeBase()`，需要同步更新构造函数调用：

```cpp
AlgoFMM::computeBase(
    octreeNodes_, wave, integralEquType_, maxLevel_,
    static_cast<int>(rwgs.size()),
    L_k1, L_k2, row, levelSpan,
    L_level_k1, L_level_k2,
    WGL_k1, WGL_k2, WGL_phi_k1, WGL_phi_k2,
    theta_level_k1, theta_level_k2,
    phi_level_k1, phi_level_k2,
    kp_lvl_k1, kp_lvl_k2);
```

FMM 单层远场仍只使用 `kp_lvl_k1[0]` / `kp_lvl_k2[0]`，所以行为主要受最细层 `L_level[0]` 影响。若希望尽量不扰动 FMM，可先让 FMM 继续使用旧采样逻辑，单独新增一个 `computeBaseMLFMM()` 给 MLFMM 使用。

推荐顺序：

1. 先给 MLFMM 新增独立的 `computeBaseMLFMM()`，降低对已验证 FMM 的影响。
2. MLFMM 验证通过后，再考虑是否统一 FMM/MLFMM 的角向采样实现。

## 11. 可选：转移核中的 `kdint` 限制

Fortran `cal_transf_kernel()` 里还有一层限制：

```fortran
Lmax=num_theta_multipole
kdint=ceiling(k_rho)
if(kdint<Lmax)then
    Lmax=kdint
end if
```

如果要完全贴近 Fortran，可以在 C++ 中加入：

```cpp
int L_eff = L;
const int kdint = static_cast<int>(std::ceil(std::abs(kR)));
if (kdint < L_eff) {
    L_eff = kdint;
}

ComplexBessel::spherical_hankel_2_array(L_eff, kR, h2_array);
for (int l = 0; l <= L_eff; ++l) {
    ...
}
```

但这一步建议作为第二阶段验证。第一阶段先把“每层 L”改正确，避免一次改变太多变量。

## 12. 验证方案

不要直接用最终 RCS 判断。建议按下面顺序验证：

1. 打印 C++ 每层 `L/theta/phi`，与 Fortran 日志对齐。

   对 `Sphere_Die_1e9`，Fortran 日志是：

   ```text
   第4层 k1: theta=16, phi=31
   第4层 k2: theta=23, phi=45
   第3层 k1: theta=23, phi=45
   第3层 k2: theta=35, phi=69
   第2层 k1: theta=35, phi=69
   第2层 k2: theta=58, phi=115
   ```

2. 对同一模型、同一随机向量 `x`，比较：

   ```text
   ||A_mlfmm*x - A_fmm*x|| / ||A_fmm*x||
   ```

3. 分别测试：

   - 只改每层 `L`。
   - 每层 `L` + Fortran 式 `theta=L+1, phi=2L+1`。
   - 再加入 `l2l` 相移符号修正。

4. 最后再跑 GMRES，看残差曲线是否从停滞转为稳定下降。

## 13. 风险点

- 角向采样点数改为 `2L+1` 后，`phi` 不再是偶数点数，所有 `idx = t * pN + p` 的代码是安全的，但任何隐含偶数假设都需要检查。
- `computeInterpolationMatrix()` 当前对 `theta/phi` 都周期回绕，这与 Fortran 移植一致；本方案先不改变它。
- `L_level` 层序必须与 `kp_lvl/theta_level/phi_level` 一致，否则转移因子会使用错误层的截断阶数。
- 如果同时改 FMM 和 MLFMM，可能影响原本已验证的 FMM 结果；建议优先只改 MLFMM 路径。

## 14. 推荐实施顺序

1. 在 `MLFMM.h` 增加 `L_level_k1/L_level_k2`。
2. 新增 `AlgoFMM::computeBaseMLFMM()`，按 Fortran 每层盒长生成 `L/theta/phi/kp`。
3. 修改 `AlgoMLFMM::computeTransferFactorMatrix()`，接收 `L_level` 并逐层使用。
4. 修改 `MLFMM.cpp` 调用点。
5. 编译并打印每层点数，先与 Fortran 日志对齐。
6. 做 `FMM Ax` 与 `MLFMM Ax` 对比。
7. 若误差仍大，再修正 `ProcessMLFMM::l2l()` 的相移符号，并继续比较。

