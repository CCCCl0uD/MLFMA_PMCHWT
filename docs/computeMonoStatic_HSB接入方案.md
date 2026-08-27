# computeMonoStatic_HSB 接入方案

## 1. 目标

当前 PMCHWT 单站 RCS 计算仍然采用逐角度直接求解：

```cpp
RCSUtils::computeMonoStatic_PMCHWT(*this, cfg, computeV);
```

该方式会对扫描角域中的每一个角度重新计算右端项、重新求解一次线性方程组。当扫描点很多，尤其是复杂介质目标或高对比度介质目标时，计算量会很大。

`computeMonoStatic_HSB` 的目标是引入空间带宽有限性，先根据目标电尺寸建立有限角域采样库，只在采样方向上求解，然后对目标输出角度进行插值，从而减少单站宽角域扫描中的重复求解次数。

## 2. 理论依据

根据空间带宽有限性，目标远区散射场可用有限阶球谐或角域有限带宽函数近似表示。设目标被半径为 `a` 的球包围，工作波数为

```text
k = 2*pi/lambda
```

则空间带宽近似为：

```text
W = rho * k * a
```

其中 `rho` 是超带宽因子，一般取大于 1 的安全系数。当前代码默认：

```cpp
rho = 1.2
```

俯仰方向采样数：

```text
M = ceil(2W)
```

第 `m` 个俯仰采样角：

```text
theta_m = m*pi/(M+1),  m = 1, 2, ..., M
```

给定 `theta_m` 时，方位方向有效带宽：

```text
N_m = ceil(W * sin(theta_m))
```

方位采样角：

```text
phi_mn = 2*pi*n/(2N_m + 1),  n = -N_m, ..., N_m
```

在这些有限采样点上直接调用 MLFMA/PMCHWT 求解，得到复散射场向量。对非采样角度，通过角域插值得到复散射场，再投影成垂直/水平极化 RCS。

## 3. 当前代码状态

### 3.1 已经具备的函数

HSB 单站函数已经存在于：

```text
Code/RCS.h
```

函数入口：

```cpp
template<typename SolverType, typename ComputeVFunc>
inline void computeMonoStatic_HSB(SolverType& solver, const RCSExportConfig& cfg,
    ComputeVFunc computeV, double rho = 1.2)
```

它目前已经完成以下工作：

1. 用 `estimateBoundingSphereRadius(solver)` 估计目标包围球半径 `a`。
2. 用 `W = rho * k1 * radius` 计算空间带宽。
3. 按 `theta_m` 和 `phi_mn` 建立二维角域采样点。
4. 在每个采样方向上调用 `computeV` 和 `solver.matrix_solver`。
5. 对 PMCHWT 和 PEC 情况分别调用 `calculateRCS_PMCHWT` 或 `calculateRCS`。
6. 对用户指定的输出角度调用 `interpolateHSBField` 插值。
7. 输出 `_RCS_HSB.txt`。

### 3.2 当前没有接入的地方

PMCHWT 单站入口在：

```text
Code/MLFMM.cpp
```

当前代码为：

```cpp
void MLFMM::mlfmm_Mono_Die_Pmchwt(const RCSExportConfig& cfg, const std::string pol_wave)
{
    auto computeV = [](MLFMM& solver, double kInc[3], double eInc[3], double hInc[3]) {
        RHS::computeV_PMCHWT(solver.rwgs, solver.gausspoint,
            solver.wave.k1(), kInc, eInc, hInc, solver.Vm);
        };
    RCSUtils::computeMonoStatic_PMCHWT(*this, cfg, computeV);
}
```

也就是说，虽然 `computeMonoStatic_HSB` 已经存在，但主流程没有调用它。

## 4. 最小接入方案

最小改动是直接将 PMCHWT 单站入口从普通逐角计算切换为 HSB：

```cpp
void MLFMM::mlfmm_Mono_Die_Pmchwt(const RCSExportConfig& cfg, const std::string pol_wave)
{
    auto computeV = [](MLFMM& solver, double kInc[3], double eInc[3], double hInc[3]) {
        RHS::computeV_PMCHWT(solver.rwgs, solver.gausspoint,
            solver.wave.k1(), kInc, eInc, hInc, solver.Vm);
        };
    RCSUtils::computeMonoStatic_HSB(*this, cfg, computeV);
}
```

该方案的优点是改动最少，能快速验证 HSB 方法是否可用。

缺点是普通单站和 HSB 单站不能通过配置切换，需要手动改代码。

## 5. 推荐接入方案

建议增加一个开关，使普通逐角求解和 HSB 求解都保留。

### 5.1 在配置结构中增加开关

建议在 `RCSExportConfig` 中增加：

```cpp
bool useHSBMono = false;
double hsbRho = 1.2;
```

含义：

```text
useHSBMono = false  使用原始逐角单站求解
useHSBMono = true   使用空间带宽有限性 HSB 单站求解
hsbRho              空间带宽超带宽因子
```

### 5.2 在 batch_config 中增加字段

建议在批处理配置中支持：

```text
use_hsb_mono = 1
hsb_rho = 1.2
```

解析时映射到：

```cpp
cfg.useHSBMono = true;
cfg.hsbRho = 1.2;
```

### 5.3 修改 PMCHWT 单站入口

推荐写法：

```cpp
void MLFMM::mlfmm_Mono_Die_Pmchwt(const RCSExportConfig& cfg, const std::string pol_wave)
{
    auto computeV = [](MLFMM& solver, double kInc[3], double eInc[3], double hInc[3]) {
        RHS::computeV_PMCHWT(solver.rwgs, solver.gausspoint,
            solver.wave.k1(), kInc, eInc, hInc, solver.Vm);
        };

    if (cfg.useHSBMono) {
        RCSUtils::computeMonoStatic_HSB(*this, cfg, computeV, cfg.hsbRho);
    }
    else {
        RCSUtils::computeMonoStatic_PMCHWT(*this, cfg, computeV);
    }
}
```

如果后续希望 PEC 单站也使用 HSB，可以在 `mlfmm_Mono_Pec_Efie` 和 `mlfmm_Mono_Pec_Cfie` 中使用同样的分支。

## 6. 必须注意的问题

### 6.1 HSB 当前适合宽角域二维建库

当前 `computeMonoStatic_HSB` 是二维球面建库：

```text
theta_m: 1 ... M
phi_mn: -N_m ... N_m
```

如果用户只扫描一个固定平面，例如：

```text
theta = 90 deg, phi = 0 ... 360 deg
```

二维建库可能比逐角求解还贵。此时应考虑一维 HSB 版本，只沿扫描角度建库和插值。

建议判断：

```text
如果 scaThetaStart == scaThetaEnd，则是固定 theta 的 phi 扫描。
如果 scaPhiStart == scaPhiEnd，则是固定 phi 的 theta 扫描。
```

固定平面扫描时，推荐后续增加：

```cpp
computeMonoStatic_HSB_1D(...)
```

### 6.2 当前 theta 插值核需要进一步校准

当前插值函数位于：

```text
Code/OverloadAlgo.h
```

核心函数：

```cpp
dirichletKernel(...)
interpolateHSBField(...)
```

当前实现采用周期 Dirichlet 核：

```cpp
sin((2*n + 1)*delta/2) / ((2*n + 1)*sin(delta/2))
```

它对方位方向比较自然，因为 `phi_mn` 是周期均匀采样。但是对俯仰方向，采样点是：

```text
theta_m = m*pi/(M+1)
```

不是 `0 ... 2pi` 上的周期均匀采样。因此当前 theta 插值不一定严格满足采样点处精确重构。初步验证可以先使用，但正式使用前建议改成更匹配该节点的 sinc/正弦级数插值。

### 6.3 应插值复场，不应插值 dB RCS

当前代码插值的是 `RCS_pre` 复远场三分量：

```cpp
sample.field = std::move(RCS_pre);
```

这是正确的。不要先计算 dB RCS 再插值，否则相位信息丢失，峰谷位置容易失真。

### 6.4 PMCHWT 情况下应复用新预条件器

`computeMonoStatic_HSB` 内部调用：

```cpp
solver.matrix_solver(...)
```

因此如果 `MLFMM::matrix_solver` 已经为 PMCHWT 接入 `2x2 block Jacobi` 或其他预条件器，HSB 建库求解会自动受益。

## 7. 建议验证流程

### 7.1 输出直接解作为基准

先保持：

```cpp
RCSUtils::computeMonoStatic_PMCHWT(*this, cfg, computeV);
```

得到普通直接求解结果：

```text
*_RCS.txt
```

### 7.2 启用 HSB 输出

切换为：

```cpp
RCSUtils::computeMonoStatic_HSB(*this, cfg, computeV, 1.2);
```

得到：

```text
*_RCS_HSB.txt
```

### 7.3 对比误差

对相同扫描角度比较：

```text
RCS_direct(theta, phi)
RCS_HSB(theta, phi)
```

建议统计：

```text
最大绝对误差: max |RCS_HSB - RCS_direct|
均方根误差:   sqrt(mean((RCS_HSB - RCS_direct)^2))
峰值位置误差: 主瓣/尖峰角度是否偏移
```

可接受误差需根据论文或工程需求确定。初期建议目标：

```text
大部分角度误差 < 1 dB
主瓣峰值误差 < 0.5 dB
深零点附近允许误差较大
```

### 7.4 扫描 rho 敏感性

建议测试：

```text
rho = 1.0
rho = 1.2
rho = 1.5
rho = 2.0
```

记录：

```text
HSB sample points
总求解时间
相对直接解误差
```

如果 `rho` 增大后误差明显下降，说明原采样带宽不足。如果 `rho` 增大后误差变化很小，则误差主要可能来自插值核或数值求解误差。

## 8. 建议实现顺序

### 阶段一：最小接入

直接将 PMCHWT 单站入口改为：

```cpp
RCSUtils::computeMonoStatic_HSB(*this, cfg, computeV);
```

完成编译，跑一个小目标算例确认能输出 `_RCS_HSB.txt`。

### 阶段二：加入配置开关

增加：

```cpp
cfg.useHSBMono
cfg.hsbRho
```

并在 `mlfmm_Mono_Die_Pmchwt` 中分支调用。

### 阶段三：误差评估

保留直接求解结果和 HSB 结果，写 MATLAB 或 Python 脚本对比 RCS 曲线。

### 阶段四：优化插值核和一维扫描

若误差不满足要求：

1. 先提高 `rho`。
2. 再修正 theta 插值核。
3. 如果只做固定平面扫描，实现一维 HSB 插值，减少建库点数。

## 9. 最终推荐方案

推荐不要永久替换 `computeMonoStatic_PMCHWT`，而是保留两条路径：

```text
普通逐角求解：精度基准、调试、论文对比
HSB 快速求解：宽角域扫描、批量计算、运动目标角域库
```

最终入口建议为：

```cpp
if (cfg.useHSBMono) {
    RCSUtils::computeMonoStatic_HSB(*this, cfg, computeV, cfg.hsbRho);
}
else {
    RCSUtils::computeMonoStatic_PMCHWT(*this, cfg, computeV);
}
```

这样可以在不破坏原有 PMCHWT 单站结果的前提下，引入空间带宽有限性，并方便后续做精度和加速比对比。
