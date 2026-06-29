# PMCHWT MoM 不收敛问题排查记录

## 结论摘要

按你的要求，这份记录只关注稠密 MoM，不讨论 FMM 和 MLFMM。

当前 PMCHWT 在 MoM 下 GMRES 残差停在约 `0.73`、CGS 不收敛，最可疑的位置不是迭代器，而是 MoM 的 PMCHWT 系统矩阵和右端向量没有使用一致的未知量定义、阻抗尺度和块矩阵符号。

最高优先级怀疑点：

1. `Code/Zmartix.h::computeZ_die_pmchwt` 的四块 PMCHWT 装配公式与参考 MoM 工程不一致。
2. `Code/Vm.h::computeV_PMCHWT` 的第二块右端没有乘阻抗尺度，而参考工程使用 `eta0 * H_inc`。
3. 当前 `OLK::L_operator` 是无阻抗 L 算子，但 PMCHWT 装配时混用了 `eta1`、`eta2`、`1/eta1`、`1/eta2`，这套归一化需要重新核对。
4. 介质损耗虚部符号可能与 `exp(-j k R)` 格林函数约定不匹配，这会增加病态性，但优先级低于矩阵块公式。

## 当前 MoM 代码路径

位置：`Code/MoM.cpp`，约第 45-47 行。

MoM 的 PMCHWT 路径是：

```cpp
if (integralEquType_ == 2) {
    Zmartix::computeZ_die_pmchwt(rwgs, row, wave, gausspoint, Z_mom);
}
```

迭代器只是在解：

```text
Z_mom * x = Vm
```

因此 MoM 下的关键只剩两个东西：

- `Z_mom` 是否正确表达 PMCHWT 边界条件。
- `Vm` 是否与 `Z_mom` 的未知量尺度一致。

## 当前 MoM PMCHWT 矩阵装配

位置：`Code/Zmartix.h::computeZ_die_pmchwt`，约第 511-520 行。

当前代码先计算无阻抗版 L/K 算子：

```cpp
OLK::L_operator(rwgField, rwgSource, wave.k1(), gausspoint, ZL1);
OLK::L_operator(rwgField, rwgSource, wave.k2(), gausspoint, ZL2);

OLK::K_operator(rwgField, rwgSource, wave.k1(), gausspoint, ZK1);
OLK::K_operator(rwgField, rwgSource, wave.k2(), gausspoint, ZK2);
```

然后装配：

```cpp
z11 = eta1 * ZL1 - eta2 * ZL2;
z12 = ZK1 - ZK2;
z21 = -ZK1 + ZK2;
z22 = (1 / eta1) * ZL1 - (1 / eta2) * ZL2;
```

这个形式很可疑：

- 内外区域贡献是相减：`L1 - L2`、`K1 - K2`。
- `z12/z21` 没有阻抗尺度。
- `z22` 使用 `1/eta` 缩放，和参考 MoM 工程的 `ZL1 + epsilon_r * ZL2` 型式不一致。
- 如果未知量第二块是磁流 `M`，那么四块矩阵必须明确对应同一种 `M` 的定义，否则即使单独 L/K 都正确，组合后的 PMCHWT 也会错。

## 与 Integral_Solver_V3.2_OMP 的 MoM 对照

参考位置：

`D:\MyCode\Integral_Solver_V3.2_OMP\Integral_Solver_V3.2_OMP\Integral Solver Kernel\MoM\MoM.f90`

`assemble_Z_PMCHW` 约第 1131-1135 行：

```fortran
this%Cz(edge_f_ID, edge_s_ID)            = ZL1 + ZL2
this%Cz(edge_f_ID, edge_s_ID + RWGNum)   = etar1 * (ZK1 + ZK2)
this%Cz(edge_f_ID + RWGNum, edge_s_ID)   = -etar1 * (ZK1 + ZK2)
this%Cz(edge_f_ID + RWGNum, edge_s_ID + RWGNum) = ZL1 + ZL2 * csr
```

其中：

```fortran
ZL1 = cal_EFIE_kernel(edge_f_ID, edge_s_ID, Obj, k1, etar1)
ZL2 = cal_EFIE_kernel(edge_f_ID, edge_s_ID, Obj, k2, etar2)
ZK1 = cal_MFIE_kernel(edge_f_ID, edge_s_ID, Obj, k1, .true.)
ZK2 = cal_MFIE_kernel(edge_f_ID, edge_s_ID, Obj, k2, .true.)
```

也就是说，参考 MoM 工程是：

```text
Z11 = ZL1 + ZL2
Z12 = eta1 * (ZK1 + ZK2)
Z21 = -eta1 * (ZK1 + ZK2)
Z22 = ZL1 + epsilon_r * ZL2
```

这与当前工程的：

```text
Z11 = eta1 * L1 - eta2 * L2
Z12 = K1 - K2
Z21 = -K1 + K2
Z22 = L1/eta1 - L2/eta2
```

差异非常大。这个差异足够导致 GMRES 残差停在较高水平。

## OLK::L_operator 与 EFIE::computeEFIE_Zij 的尺度关系

当前 MoM PMCHWT 使用 `OLK::L_operator`，不是 `EFIE::computeEFIE_Zij`。

从代码形式看：

`OLK::L_operator` 的核心类似：

```text
j*k * (G * dot - G/k^2 * divdiv)
```

`EFIE::computeEFIE_Zij(k, eta)` 的核心类似：

```text
j*k*eta * G * dot - j*eta/k * G * divdiv
```

因此近似关系是：

```text
EFIE_Z(k, eta) = eta * L_operator(k)
```

如果你想沿用参考工程的 MoM 块公式，同时继续使用 `OLK::L_operator`，更一致的写法应先构造已含阻抗的 EFIE 块：

```cpp
E1 = wave.eta1() * ZL1;
E2 = wave.eta2() * ZL2;
```

然后再装配：

```cpp
Z11 = E1 + E2;
Z12 = wave.eta1() * (ZK1 + ZK2);
Z21 = -wave.eta1() * (ZK1 + ZK2);
Z22 = E1 + wave.epsilonR() * E2;
```

这并不是说它已经百分百正确，因为还要看你的 `K_operator` 与参考工程 `cal_MFIE_kernel(..., .true.)` 的法向和奇异项符号是否一致。但它至少比当前的相减形式更接近参考 MoM 工程。

## 右端向量问题

位置：`Code/Vm.h::computeV_PMCHWT`，约第 181-200 行。

当前右端：

```cpp
Vm[rwgid]     += <f, E_inc>;
Vm[N + rwgid] += <f, H_inc>;
```

参考工程右端：

`D:\MyCode\Integral_Solver_V3.2_OMP\Integral_Solver_V3.2_OMP\Exitation\Excitation.f90`，约第 177-179 行：

```fortran
V(i_RWG)          = Ei
V(i_RWG + RWGNum) = etar0 * Hi
```

当前 `EMSource::initPW` 中，`vH` 是由：

```cpp
crossPW(data_vW, data_vE, data_vH);
normalize(data_vH);
```

得到的方向向量，它没有包含物理磁场幅度 `1/eta0`。因此当前 `Vm[N + rwgid]` 的尺度很可能和矩阵第二行不匹配。

如果采用参考工程的 PMCHWT 归一化，应把第二块右端改为：

```cpp
Vm[N + rwgid] += wave.eta1() * lm * gp.weight[i] * dot(hInc, rho) * e_k1;
```

正负面两处都要同步改。

## K 算子法向与介质 PMCHWT 的风险

你提到 K 算子是从 MFIE 中提取出来的。这里需要非常谨慎：PEC 的 MFIE 中经常内含 `n x`、主值项、跳跃项或特定法向方向约定；PMCHWT 的 K 块对内外区域的法向、符号、是否包含跳跃项十分敏感。

建议做一个单元级检查：

1. 选一个很小闭合网格。
2. 对同一对 RWG，分别打印 `ZK1`、`ZK2`。
3. 与 `Integral_Solver_V3.2_OMP` 的 `cal_MFIE_kernel(..., .true.)` 对比符号趋势。
4. 特别注意当前 `K_operator` 是否已经包含 `1/2 I` 跳跃项。如果 PMCHWT 块公式默认 K 是主值算子，而你的 K 已含跳跃项，矩阵会错。

CGS 对这类非正规复矩阵很脆弱，K 块符号一错通常会比 GMRES 更早表现为发散或 breakdown。

## 介质参数虚部符号

位置：`Code/main.cpp`，约第 97 行：

```cpp
std::complex<double> epsilonR(4.0, 0.001);
```

位置：`Code/EMSource.cpp::initDielectric`：

```cpp
// 有耗介质 （εr = ε' - jε"）
```

当前格林函数是 `exp(-J * k * R)`。如果代码约定有耗介质是 `epsilon = epsilon' - j epsilon''`，那么示例中的 `+0.001j` 可能对应增益而不是损耗。

建议打印：

```cpp
wave.k2()
```

并确认在 `exp(-j*k2*R)` 下场随距离衰减。如果虚部符号反了，矩阵会更病态。但我仍认为它不是第一根因，第一根因更像 PMCHWT 块装配。

## 建议验证顺序

1. 固定 `selectAlgorithm = 0`，只跑 MoM。
2. 暂时只改 `computeZ_die_pmchwt` 和 `computeV_PMCHWT`，不要碰 FMM/MLFMM。
3. 先用最小介质球网格验证 GMRES 是否能明显低于 `0.73`。
4. 比较 `epsilonR = 4.0 - 0.001j` 与 `4.0 + 0.001j` 的残差趋势，确认损耗符号。
5. 如果矩阵改完后 GMRES 收敛但 CGS 仍差，优先接受 GMRES 结果；CGS 对 PMCHWT 本来就不如 GMRES 稳。

## 最可能根因排序

1. **`computeZ_die_pmchwt` 的四块 PMCHWT 公式不符合当前未知量定义**。
2. **`computeV_PMCHWT` 第二块右端缺少阻抗尺度**。
3. **从 MFIE 提取的 `K_operator` 可能包含了不适合 PMCHWT 的法向/跳跃项约定**。
4. **介质损耗虚部符号可能反了，导致病态性变强**。
5. **CGS 算法本身不适合当前非正规复矩阵**。

我的建议是先不要调 `itr_max`、`tol` 或换迭代器。先让 MoM 的 PMCHWT 矩阵和右端在量纲、符号、未知量定义上自洽；否则 GMRES/CGS 都只是在努力求解一个错误或严重失衡的线性系统。

1.请给我这两项目使用PMCHWT公式填充阻抗矩阵的公式(包括使用算子的公式)D:\MyCode\Integral_Solver_V3.2_OMP和D:\MyCode\for90-mom2-master

2.我的computeV_PMCHWT第二块右端使用的hinc是通过电场单位向量旋转得到的
