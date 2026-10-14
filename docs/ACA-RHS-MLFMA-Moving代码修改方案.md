# ACA-RHS-MLFMA-Moving 代码修改方案

## 1. 目标

本方案用于在当前 `PMCHWT_MLFMA/MLFMA` 工程中新增一条运动 PEC 目标加速路径：

```text
ACA-RHS-MLFMA-Moving
= 运动目标角度等效 + 右边向量 ACA 压缩 + 固定 MLFMA 阻抗算子 + 逐时刻 RCS 直接计算
```

第一版的核心要求：

1. 参考现有 `HSB-MLFMA-Moving` 的运动目标处理方式。
2. 利用远场目标平移、旋转不变性，将目标运动等效为右边向量变化。
3. 保持阻抗矩阵/MLFMA 算子不变，不随时间步重建。
4. 不使用 HSB 插值，不插值 RCS。
5. 每个运动时刻仍然输出由电流合成后直接计算得到的 RCS 值。

数学形式：

```text
Z * I(t) = b(t)
```

其中：

```text
Z    固定，由 MLFMA 隐式矩阵向量乘提供
b(t) 随运动时刻的目标本体系入射角 theta_body(t), phi_body(t) 改变
I(t) 由 ACA 模态解线性合成
RCS(t) 由 I(t) 逐时刻直接调用 calculateRCS 计算
```

---

## 2. 第一版支持范围

第一版建议固定为：

```text
目标类型：PEC
算法：MLFMA
收发模式：Mono
积分方程：CFIE 优先，EFIE 可同步支持
路径输入：moving_path_file
路径格式：time_s x y z qw qx qy qz
雷达：单站 radar_position
RCS 输出：运动时间序列 RCS
```

第一版不支持：

```text
介质 PMCHWT
双站运动目标
多目标相对运动
距离相位 exp(-j*2*k*range)
HSB 插值
RCS/dB 插值
运动过程中重建几何、OCTree、Z_near、MLFMA 转移因子
```

说明：当前工程中 `IntegralEquation::PMCHWT` 对应介质目标路径。用户需求是运动 PEC，因此第一版应走 `IntegralEquation::EFIE` 或 `IntegralEquation::CFIE`，不走 PMCHWT。

---

## 3. 运动等效关系

沿用 `HSB-MLFMA-Moving` 的坐标处理。

路径文件每一行：

```text
time_s  x  y  z  qw  qx  qy  qz
```

含义：

```text
time_s        当前时刻
x y z         目标质心在全局坐标系的位置
qw qx qy qz   目标姿态四元数，表示 body -> global 的旋转
```

单站雷达位置来自配置：

```text
radar_position = radarX, radarY, radarZ
```

每个时刻的等效入射/观测方向：

```text
r_global(t) = normalize(radar_position - target_position(t))
R(t)         = body_to_global rotation matrix
r_body(t)   = R(t)^T * r_global(t)
```

再转换为目标本体系角度：

```text
theta_body(t) = acos(r_body.z)
phi_body(t)   = atan2(r_body.y, r_body.x)
```

极点规则沿用现有代码：

```text
theta = 0 或 180 时，phi = 0
```

复用现有文件：

```text
Code/MovingTarget.h
Code/MovingTarget.cpp
```

复用现有函数：

```cpp
MovingTarget::readMovingPathFile(...)
MovingTarget::computeMovingAngleSamples(...)
MovingTarget::computeMovingAngleSample(...)
MovingTarget::quaternionToMatrix(...)
MovingTarget::rotateGlobalToBody(...)
MovingTarget::vectorToThetaPhiDeg(...)
```

---

## 4. ACA 压缩对象

把所有运动时刻的右边向量组成矩阵：

```text
B in C^(N x T)
```

其中：

```text
N = RWG 个数
T = moving_path_file 中的运动时刻数
B(i,t) = 第 t 个时刻、第 i 条 RWG 的 RHS 值
```

对于 PEC EFIE：

```text
B(i,t) = V_E(i, theta_body(t), phi_body(t))
```

对于 PEC CFIE：

```text
B(i,t) = alpha * V_E(i,t) + (1-alpha) * V_M(i,t)
```

ACA 分解：

```text
B ~= U_rhs * V_rhs^T
```

存储约定建议使用 rank-major：

```text
U_rhs[m * N + i] = U(i,m)
V_rhs[m * T + t] = V(t,m)
```

线性求解：

```text
Z * y_m = U_rhs(:,m),  m = 0..rank-1
```

逐时刻电流合成：

```text
I(t) = sum_m y_m * V_rhs(m,t)
```

逐时刻 RCS 直接计算：

```text
RCS(t) = calculateRCS(I(t), theta_body(t), phi_body(t))
```

---

## 5. 新增配置项

修改文件：

```text
Code/RCSExportConfig.h
Code/BatchRunner.h
Code/main.cpp
Code/batch_config.txt
```

### 5.1 RCSExportConfig.h

在 `struct RCSExportConfig` 中新增：

```cpp
bool useACARHSMoving{ false };
double acaRHSTol{ 1.0e-3 };
int acaRHSMaxRank{ 0 };
double acaRHSPivotEps{ 1.0e-14 };
bool exportACARHSCheck{ true };
int acaRHSCheckCount{ 10 };
bool exportACARHSFactors{ false };
```

字段含义：

```text
useACARHSMoving      启用 ACA-RHS-MLFMA-Moving
acaRHSTol            ACA RHS 相对截断阈值
acaRHSMaxRank        ACA 最大秩，0 表示 min(N,T)
acaRHSPivotEps       pivot 过小时停止，防止除以近零数
exportACARHSCheck    输出 RHS 重构抽样误差
acaRHSCheckCount     抽样检查的时间步数量
exportACARHSFactors  可选输出 U/V 因子，默认关闭
```

### 5.2 validate() 规则

在 `RCSExportConfig::isValid()` 中新增：

```cpp
if (useACARHSMoving) {
    if (selectAlgo != selectAlgorithm::MLFMA) return false;
    if (material != materialType::PEC) return false;
    if (monoDual != MonoDual::Mono) return false;
    if (integralEqu == IntegralEquation::PMCHWT) return false;
    if (movingPathFile.empty()) return false;
    if (acaRHSTol <= 0.0) return false;
    if (acaRHSMaxRank < 0) return false;
    if (acaRHSPivotEps <= 0.0) return false;
    if (acaRHSCheckCount < 0) return false;
}
```

建议 `HSB-MLFMA-Moving` 和 `ACA-RHS-MLFMA-Moving` 互斥：

```cpp
if (useACARHSMoving && useHSBMLFMAMoving) return false;
```

### 5.3 输出文件名标识

在 `RCSExportConfig::toCoreString()` 中追加：

```cpp
if (useACARHSMoving) {
    oss << "_ACARHSMoving";
}
```

### 5.4 BatchRunner.h

在 `SimTask` 中新增同名字段：

```cpp
bool useACARHSMoving;
double acaRHSTol;
int acaRHSMaxRank;
double acaRHSPivotEps;
bool exportACARHSCheck;
int acaRHSCheckCount;
bool exportACARHSFactors;
```

在 `BatchConfig` 中新增默认值：

```cpp
bool useACARHSMoving = false;
double acaRHSTol = 1.0e-3;
int acaRHSMaxRank = 0;
double acaRHSPivotEps = 1.0e-14;
bool exportACARHSCheck = true;
int acaRHSCheckCount = 10;
bool exportACARHSFactors = false;
```

在 `BatchConfig::generateTasks()` 中复制字段到 `task`。

在 `BatchConfigParser::parse()` 中新增 key：

```cpp
else if (key == "use_aca_rhs_moving") {
    cfg.useACARHSMoving = parseBool(value);
}
else if (key == "aca_rhs_tol") {
    cfg.acaRHSTol = std::stod(value);
}
else if (key == "aca_rhs_max_rank") {
    cfg.acaRHSMaxRank = std::stoi(value);
}
else if (key == "aca_rhs_pivot_eps") {
    cfg.acaRHSPivotEps = std::stod(value);
}
else if (key == "export_aca_rhs_check") {
    cfg.exportACARHSCheck = parseBool(value);
}
else if (key == "aca_rhs_check_count") {
    cfg.acaRHSCheckCount = std::stoi(value);
}
else if (key == "export_aca_rhs_factors") {
    cfg.exportACARHSFactors = parseBool(value);
}
```

在 `executeSingleTask()` 中把 task 字段写入 `cfg`。

### 5.5 main.cpp

单次运行调试区新增：

```cpp
bool useACARHSMoving = false;
double acaRHSTol = 1.0e-3;
int acaRHSMaxRank = 0;
double acaRHSPivotEps = 1.0e-14;
bool exportACARHSCheck = true;
int acaRHSCheckCount = 10;
bool exportACARHSFactors = false;
```

并赋值：

```cpp
cfg.useACARHSMoving = useACARHSMoving;
cfg.acaRHSTol = acaRHSTol;
cfg.acaRHSMaxRank = acaRHSMaxRank;
cfg.acaRHSPivotEps = acaRHSPivotEps;
cfg.exportACARHSCheck = exportACARHSCheck;
cfg.acaRHSCheckCount = acaRHSCheckCount;
cfg.exportACARHSFactors = exportACARHSFactors;
```

### 5.6 batch_config.txt 示例

新增示例配置：

```text
# ACA RHS moving acceleration: 0 = off, 1 = on
use_aca_rhs_moving = 0
aca_rhs_tol = 1.0e-3
aca_rhs_max_rank = 0
aca_rhs_pivot_eps = 1.0e-14
export_aca_rhs_check = 1
aca_rhs_check_count = 10
export_aca_rhs_factors = 0
```

启用时推荐：

```text
algorithm = 2
integral_equation = 1
mono_dual = mono
use_aca_rhs_moving = 1
use_hsb_mlfma_moving = 0
moving_path_file = D:\MyCode\PMCHWT_MLFMA\DATA\moving_path.txt
radar_position = 0,0,0
```

---

## 6. 新增 ACA_RHS 模块

新增文件：

```text
Code/ACARHS.h
```

因为 ACA 模块需要接受 lambda，建议第一版做成纯头文件模板。

### 6.1 头文件结构

```cpp
#pragma once
#ifndef ACA_RHS_H
#define ACA_RHS_H

#include <algorithm>
#include <cmath>
#include <complex>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <vector>

namespace ACARHS {

struct Options {
    double tol = 1.0e-3;
    int maxRank = 0;
    double pivotEps = 1.0e-14;
};

struct Factors {
    int rows = 0;
    int cols = 0;
    int rank = 0;
    std::vector<std::complex<double>> U;
    std::vector<std::complex<double>> V;
};

inline double abs2(const std::complex<double>& z);
inline bool containsIndex(const std::vector<int>& indices, int value);

template<typename ComputeColumnFunc, typename ComputeEntryFunc>
Factors factorize(
    int rows,
    int cols,
    const Options& opt,
    ComputeColumnFunc computeColumn,
    ComputeEntryFunc computeEntry);

inline std::vector<std::complex<double>> reconstructColumn(
    const Factors& factors,
    int col);

} // namespace ACARHS

#endif
```

### 6.2 computeColumn 与 computeEntry 约定

`computeColumn`：

```cpp
computeColumn(int timeIndex, std::vector<std::complex<double>>& rhs)
```

要求：

```text
输入 timeIndex
输出完整 RHS 向量 rhs，长度 N
```

`computeEntry`：

```cpp
std::complex<double> computeEntry(int rowIndex, int timeIndex)
```

要求：

```text
只计算 B(rowIndex, timeIndex)
```

### 6.3 ACA 算法流程

伪代码：

```text
rows = N
cols = T
rankLimit = maxRank > 0 ? min(maxRank, rows, cols) : min(rows, cols)

I_k = 0
selectedRows = {0}
selectedCols = {}

while rank < rankLimit:
    1. 固定行 I_k，扫描所有时间列：
       rowResidual[t] = B(I_k,t) - sum_p U[p,I_k] * V[p,t]

    2. 选择未选过且 |rowResidual[t]| 最大的 J_k。
       若 pivot 太小，停止。

    3. 新 V 模态：
       V[rank,t] = rowResidual[t] / rowResidual[J_k]

    4. 固定列 J_k，生成完整 RHS：
       column = B(:,J_k)
       colResidual[i] = column[i] - sum_p U[p,i] * V[p,J_k]

    5. 新 U 模态：
       U[rank,i] = colResidual[i]

    6. 选择未选过且 |colResidual[i]| 最大的 I_k。

    7. 用 ||U_new||^2 * ||V_new||^2 判断相对贡献：
       如果 newContribution <= tol^2 * approxNorm，则停止。
```

第一版可以采用较稳妥的停止准则：

```cpp
const double newContribution = norm2U * norm2V;
approxNorm += newContribution;
if (rank > 0 && newContribution <= opt.tol * opt.tol * approxNorm) {
    break;
}
```

后续如需更接近旧 ACA 项目，可再补交叉项 Frobenius 更新。

---

## 7. 修改 RHS 计算接口

修改文件：

```text
Code/Vm.h
```

当前已有整列 RHS：

```cpp
RHS::computeV_E(...)
RHS::computeV_M(...)
RHS::computeV_C(...)
```

新增单个 RWG 版本：

```cpp
inline std::complex<double> computeV_E_entry(
    const std::vector<RWGBase>& rwgs,
    const gaussPoints& gp,
    const std::complex<double> k_,
    const double kInc[3],
    const double eInc[3],
    int rwgid);

inline std::complex<double> computeV_M_entry(
    const std::vector<RWGBase>& rwgs,
    const gaussPoints& gp,
    const std::complex<double> k_,
    const double kInc[3],
    const double hInc[3],
    int rwgid);

inline std::complex<double> computeV_C_entry(
    const std::vector<RWGBase>& rwgs,
    const gaussPoints& gp,
    const std::complex<double> k_,
    const double kInc[3],
    const double eInc[3],
    const double hInc[3],
    const double alpha_,
    int rwgid);
```

实现要求：

1. 从现有 `computeV_E` 和 `computeV_M` 中抽取单个 `rwgid` 的高斯积分逻辑。
2. `computeV_E` 可以继续保留原实现；也可以内部循环调用 `computeV_E_entry`，但需注意性能。
3. `computeV_C_entry` 返回：

```cpp
alpha_ * computeV_E_entry(...) + (1.0 - alpha_) * computeV_M_entry(...)
```

注意：

```text
ACA-RHS-Moving 第一版只用于 PEC EFIE/CFIE，因此不需要新增 PMCHWT entry。
```

---

## 8. 新增 RCSUtils 主流程

修改文件：

```text
Code/RCS.h
```

在现有 `computeMonoStatic_HSB_Moving` 附近新增：

```cpp
template<typename SolverType, typename ComputeVFunc, typename ComputeVEntryFunc>
inline void computeMonoStatic_ACARHS_Moving(
    SolverType& solver,
    const RCSExportConfig& cfg,
    ComputeVFunc computeV,
    ComputeVEntryFunc computeVEntry)
```

需要包含：

```cpp
#include "ACARHS.h"
```

### 8.1 函数流程

```text
1. 检查 solver.integralEquType_ != 2。
2. 读取 cfg.movingPathFile。
3. 使用 MovingTarget 计算 movingAngles。
4. 可复用或新增导出角度映射文件。
5. 构造 computeColumn(timeIndex, rhs)。
6. 构造 computeEntry(rowIndex, timeIndex)。
7. 调用 ACARHS::factorize 得到 U/V。
8. 对每个 rank 模态调用 solver.matrix_solver。
9. 对每个原始运动时刻合成 I(t)。
10. 对每个时刻直接 calculateRCS。
11. projectRCSGridValue 得到 rcs_db / rcs_v_db / rcs_h_db。
12. 输出 *_RCS_ACARHS_MLFMA_MOVING_*.txt。
13. 可选输出 ACA RHS 重构误差检查。
14. 可选输出 ACA 因子。
```

### 8.2 入射角生成

对第 `t` 个运动时刻：

```cpp
const auto& angle = movingAngles[t];
EMSource pw;
pw.initPW(data_pol, angle.thetaBodyDeg, angle.phiBodyDeg);

double kinc[3] = { pw.vW(0), pw.vW(1), pw.vW(2) };
double einc[3] = { pw.vE(0), pw.vE(1), pw.vE(2) };
double hinc[3] = { pw.vH(0), pw.vH(1), pw.vH(2) };
```

### 8.3 computeColumn

```cpp
auto computeColumn = [&](int timeIndex,
    std::vector<std::complex<double>>& rhs) {
    const auto& angle = movingAngles[timeIndex];
    EMSource pw;
    pw.initPW(data_pol, angle.thetaBodyDeg, angle.phiBodyDeg);
    double kinc[3] = { pw.vW(0), pw.vW(1), pw.vW(2) };
    double einc[3] = { pw.vE(0), pw.vE(1), pw.vE(2) };
    double hinc[3] = { pw.vH(0), pw.vH(1), pw.vH(2) };
    computeV(solver, kinc, einc, hinc);
    rhs.assign(solver.Vm.begin(), solver.Vm.end());
};
```

### 8.4 computeEntry

```cpp
auto computeEntry = [&](int rowIndex, int timeIndex) {
    const auto& angle = movingAngles[timeIndex];
    EMSource pw;
    pw.initPW(data_pol, angle.thetaBodyDeg, angle.phiBodyDeg);
    double kinc[3] = { pw.vW(0), pw.vW(1), pw.vW(2) };
    double einc[3] = { pw.vE(0), pw.vE(1), pw.vE(2) };
    double hinc[3] = { pw.vH(0), pw.vH(1), pw.vH(2) };
    return computeVEntry(solver, kinc, einc, hinc, rowIndex);
};
```

### 8.5 解 ACA 模态

```cpp
std::vector<std::complex<double>> modeSolutions(
    factors.rank * solver.row,
    std::complex<double>(0.0, 0.0));

for (int m = 0; m < factors.rank; ++m) {
    std::vector<std::complex<double>> rhsMode(solver.row);
    std::vector<std::complex<double>> xMode(
        solver.row, std::complex<double>(0.0, 0.0));

    for (int i = 0; i < solver.row; ++i) {
        rhsMode[i] = factors.U[m * solver.row + i];
    }

    int itr_max = solver.row;
    int mr = 50;
    double tol_abs = 1.0e-8;
    double tol_rel = 1.0e-4;
    solver.matrix_solver(
        solver.row,
        xMode.data(),
        rhsMode.data(),
        itr_max,
        mr,
        tol_abs,
        tol_rel);

    for (int i = 0; i < solver.row; ++i) {
        modeSolutions[m * solver.row + i] = xMode[i];
    }
}
```

### 8.6 合成每个时刻的电流

```cpp
std::vector<std::complex<double>> I_t(
    solver.row, std::complex<double>(0.0, 0.0));

for (int i = 0; i < solver.row; ++i) {
    I_t[i] = std::complex<double>(0.0, 0.0);
}

for (int m = 0; m < factors.rank; ++m) {
    const std::complex<double> coeff = factors.V[m * timeCount + t];
    for (int i = 0; i < solver.row; ++i) {
        I_t[i] += modeSolutions[m * solver.row + i] * coeff;
    }
}
```

### 8.7 逐时刻直接计算 RCS

```cpp
std::vector<std::complex<double>> RCS_pre(3, 0.0);
RCSUtils::calculateRCS(
    solver,
    I_t.data(),
    RCS_pre,
    solver.wave.k1(),
    solver.wave.eta1(),
    angle.thetaBodyDeg,
    angle.phiBodyDeg);

RCSGridValue value = projectRCSGridValue(
    RCS_pre,
    data_pol,
    angle.thetaBodyDeg,
    angle.phiBodyDeg);
```

注意：

```text
这里没有 HSB 插值。
这里没有 RCS 插值。
这里对每个时刻都调用 calculateRCS。
```

---

## 9. RCS 输出设计

主输出文件：

```text
*_RCS_ACARHS_MLFMA_MOVING_*.txt
```

建议表头：

```text
idx
time_s
range_m
theta_global_deg
phi_global_deg
theta_body_deg
phi_body_deg
aca_rank
rcs_db
rcs_v_db
rcs_h_db
```

示例写法：

```cpp
rcsFile << "idx\ttime_s\trange_m\t"
    << "theta_global_deg\tphi_global_deg\t"
    << "theta_body_deg\tphi_body_deg\t"
    << "aca_rank\trcs_db\trcs_v_db\trcs_h_db\n";
```

每个时刻输出：

```cpp
rcsFile << i << "\t"
    << angle.time_s << "\t"
    << angle.range << "\t"
    << angle.thetaGlobalDeg << "\t"
    << angle.phiGlobalDeg << "\t"
    << angle.thetaBodyDeg << "\t"
    << angle.phiBodyDeg << "\t"
    << factors.rank << "\t"
    << value.rcs << "\t"
    << value.rcs_v << "\t"
    << value.rcs_h << "\n";
```

---

## 10. ACA RHS 检查输出

可选输出：

```text
*_ACA_RHS_CHECK_*.txt
```

建议表头：

```text
time_idx
time_s
rhs_exact_norm
rhs_error_norm
relative_error
```

抽样策略：

```text
如果 acaRHSCheckCount <= 0，不输出。
否则在 [0, T-1] 中均匀抽取 acaRHSCheckCount 个时间步。
```

误差计算：

```text
b_exact = computeColumn(t)
b_aca   = ACARHS::reconstructColumn(factors, t)
relative_error = ||b_exact - b_aca||_2 / ||b_exact||_2
```

注意：

```text
这个检查只验证 RHS ACA 重构误差。
RCS 误差仍需用普通逐时刻 MLFMA 基准结果比较。
```

---

## 11. ACA 因子输出

可选输出：

```text
*_ACA_RHS_FACTORS_U_*.txt
*_ACA_RHS_FACTORS_V_*.txt
```

默认关闭，因为大型算例下文件可能很大。

`U` 表头：

```text
rank_idx
row_idx
real
imag
```

`V` 表头：

```text
rank_idx
time_idx
real
imag
```

---

## 12. MLFMM.cpp 接入方案

修改文件：

```text
Code/MLFMM.cpp
```

### 12.1 PEC EFIE mono

当前位置：

```cpp
void MLFMM::mlfmm_Mono_Pec_Efie(
    const RCSExportConfig& cfg,
    const std::string pol_wave)
```

现有 `computeV`：

```cpp
auto computeV = [](MLFMM& solver, double kInc[3], double eInc[3], double hInc[3]) {
    (void)hInc;
    RHS::computeV_E(solver.rwgs, solver.gausspoint,
        solver.wave.k1(), kInc, eInc, solver.Vm);
};
```

新增 `computeVEntry`：

```cpp
auto computeVEntry = [](MLFMM& solver,
    double kInc[3], double eInc[3], double hInc[3], int rwgid) {
    (void)hInc;
    return RHS::computeV_E_entry(
        solver.rwgs,
        solver.gausspoint,
        solver.wave.k1(),
        kInc,
        eInc,
        rwgid);
};
```

分支：

```cpp
if (cfg.useACARHSMoving) {
    RCSUtils::computeMonoStatic_ACARHS_Moving(
        *this, cfg, computeV, computeVEntry);
}
else if (cfg.useHSBMLFMAMoving) {
    RCSUtils::computeMonoStatic_HSB_Moving(*this, cfg, computeV, cfg.hsbRho);
}
else {
    RCSUtils::computeMonoStatic(*this, cfg, computeV);
}
```

### 12.2 PEC CFIE mono

当前位置：

```cpp
void MLFMM::mlfmm_Mono_Pec_Cfie(
    const RCSExportConfig& cfg,
    const std::string pol_wave)
```

现有 `computeV`：

```cpp
auto computeV = [](MLFMM& solver, double kInc[3], double eInc[3], double hInc[3]) {
    RHS::computeV_C(solver.rwgs, solver.gausspoint,
        solver.wave.k1(), kInc, eInc, hInc, alpha, solver.Vm);
};
```

新增 `computeVEntry`：

```cpp
auto computeVEntry = [](MLFMM& solver,
    double kInc[3], double eInc[3], double hInc[3], int rwgid) {
    return RHS::computeV_C_entry(
        solver.rwgs,
        solver.gausspoint,
        solver.wave.k1(),
        kInc,
        eInc,
        hInc,
        alpha,
        rwgid);
};
```

推荐分支优先级：

```cpp
if (cfg.useACARHSMoving) {
    RCSUtils::computeMonoStatic_ACARHS_Moving(
        *this, cfg, computeV, computeVEntry);
}
else if (cfg.useHSBMLFMAMoving) {
    RCSUtils::computeMonoStatic_HSB_Moving(*this, cfg, computeV, cfg.hsbRho);
}
else if (cfg.useHSBMono) {
    RCSUtils::computeMonoStatic_HSB(*this, cfg, computeV, cfg.hsbRho);
}
else {
    RCSUtils::computeMonoStatic(*this, cfg, computeV);
}
```

---

## 13. Visual Studio 项目文件

如果 `ACARHS.h` 是纯头文件，需要加入：

```text
Code/MLFMA.vcxproj
Code/MLFMA.vcxproj.filters
```

在 `ClInclude` 中新增：

```xml
<ClInclude Include="ACARHS.h" />
```

如果后续拆出 `ACARHS.cpp`，再额外加入：

```xml
<ClCompile Include="ACARHS.cpp" />
```

第一版建议不建 `.cpp`。

---

## 14. 与 HSB-MLFMA-Moving 的关系

两条路径共用：

```text
moving_path_file
radar_position
MovingTarget 角度映射
角度映射检查输出逻辑
目标本体系 theta/phi 作为单站角
```

两条路径不同：

| 项目 | HSB-MLFMA-Moving | ACA-RHS-MLFMA-Moving |
|---|---|---|
| 加速对象 | 远场散射场角域采样 | 右边向量矩阵 |
| 是否插值 | 插值复远场 | 不插值 |
| 求解次数 | HSB 样本点数 | ACA rank |
| 每时刻 RCS | 由插值场投影 | 由合成电流直接 calculateRCS |
| 阻抗算子 | 固定 MLFMA | 固定 MLFMA |

优先级建议：

```text
use_aca_rhs_moving 与 use_hsb_mlfma_moving 互斥
```

原因：二者是两种不同加速路线，第一版不要叠加，以便验证误差来源。

---

## 15. 关键风险与处理

### 15.1 RHS 低秩性不足

如果运动路径覆盖角域很广，或姿态快速跳变，`B(edge,time)` 的数值秩可能升高。

处理：

```text
输出 aca_rank
输出 RHS 抽样重构误差
允许设置 aca_rhs_max_rank
与普通逐时刻 MLFMA 结果对比
```

### 15.2 pivot 接近零

ACA 归一化时如果 pivot 太小会数值不稳定。

处理：

```text
使用 aca_rhs_pivot_eps
pivot 小于阈值时停止
日志输出停止原因
```

### 15.3 RCS 误差不完全等于 RHS 误差

即使 RHS 重构误差满足阈值，RCS dB 曲线在深零点附近也可能较敏感。

处理：

```text
对比普通逐时刻 MLFMA 的 RCS
重点统计主瓣、峰值和一般角度误差
深零点单独标注
```

### 15.4 Solver 初值

每个 ACA 模态目前建议零初值求解。

后续优化：

```text
可复用前一模态解作为初值
或加入预条件器
```

第一版不做，避免引入额外变量。

---

## 16. 验证方案

### 16.1 编译验证

```text
Release x64 编译通过
use_aca_rhs_moving = 0 时旧路径结果不变
```

### 16.2 运动角度映射验证

使用单位四元数：

```text
qw qx qy qz = 1 0 0 0
```

预期：

```text
theta_global == theta_body
phi_global == phi_body
```

### 16.3 固定目标验证

路径中所有时刻位置和姿态完全相同。

预期：

```text
ACA rank 接近 1
所有时刻 RHS 几乎相同
所有时刻 RCS 相同
```

### 16.4 纯旋转验证

目标位置固定，绕 z 轴匀速旋转。

预期：

```text
theta_body 基本不变
phi_body 随姿态反向变化
RCS 随姿态周期变化
```

### 16.5 与普通逐时刻 MLFMA 对比

新增或临时保留一个 direct moving 参考流程：

```text
对每个 time step:
    compute RHS(theta_body, phi_body)
    solver.matrix_solver
    calculateRCS
```

抽取 5-10 个时刻，比较：

```text
RCS_direct(t)
RCS_ACARHS(t)
abs error in dB
```

建议初始接受标准：

```text
RHS relative error 接近 aca_rhs_tol 量级
大多数非深零点 RCS 误差 < 0.5~1 dB
主瓣峰值误差 < 0.5 dB
求解次数从 T 降为 rank
```

---

## 17. 建议实现顺序

建议按以下顺序执行，便于每一步单独编译和定位问题：

1. 修改 `RCSExportConfig.h`，加入配置字段和校验。
2. 修改 `BatchRunner.h` 和 `main.cpp`，加入配置传递。
3. 修改 `batch_config.txt`，加入示例 key。
4. 修改 `Vm.h`，新增 `computeV_E_entry/M_entry/C_entry`。
5. 新增 `ACARHS.h`，实现 ACA 分解和列重构。
6. 修改 `RCS.h`，新增 `computeMonoStatic_ACARHS_Moving`。
7. 修改 `MLFMM.cpp`，在 PEC mono EFIE/CFIE 接入新分支。
8. 修改 `MLFMA.vcxproj` 和 `.filters`，加入 `ACARHS.h`。
9. 编译。
10. 用固定目标路径验证 rank 和 RCS。
11. 用短运动路径对比普通逐时刻 MLFMA。
12. 扩展到长时间序列，统计 rank、耗时和误差。

---

## 18. 第一版最终数据流

```text
NAS 网格
  |
  v
RWG + OCTree + MLFMA 预计算
  |
  v
固定阻抗算子 Z_MLFMA
  |
  +-----------------------------+
                                |
moving_path_file                 |
  |                              |
  v                              |
MovingTarget 角度映射            |
  |                              |
  v                              |
theta_body(t), phi_body(t)       |
  |                              |
  v                              |
RHS 矩阵 B(edge,time)             |
  |                              |
  v                              |
ACA: B ~= U_rhs * V_rhs^T         |
  |                              |
  v                              |
rank 次求解 Z*y_m = U_rhs(:,m) ---+
  |
  v
每时刻合成 I(t)
  |
  v
每时刻直接 calculateRCS
  |
  v
*_RCS_ACARHS_MLFMA_MOVING_*.txt
```

---

## 19. 审核确认项

请审核以下设计选择：

1. 第一版只支持 `MLFMA + PEC + Mono + EFIE/CFIE`，不支持介质 PMCHWT。
2. 运动目标处理完全复用 `HSB-MLFMA-Moving` 的 body-frame 角度映射。
3. 阻抗矩阵/MLFMA 算子固定，不随时间步重建。
4. ACA 只压缩 RHS 矩阵 `B(edge,time)`。
5. 每个时刻 RCS 由合成电流 `I(t)` 直接调用 `calculateRCS` 得到。
6. 第一版不叠加 HSB，不插值远场，不插值 RCS。
7. 第一版不乘距离相位，只输出 RCS 计算值。
8. `use_aca_rhs_moving` 与 `use_hsb_mlfma_moving` 互斥。
9. 新增 `Code/ACARHS.h` 为纯头文件模板。
10. 新增 RHS entry 函数用于 ACA 行扫描，以避免大量完整 RHS 重算。

审核通过后，再按第 17 节顺序进行代码修改。
