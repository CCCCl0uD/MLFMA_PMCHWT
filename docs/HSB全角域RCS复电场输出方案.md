# HSB 全角域 RCS 和复电场输出修改方案

## 1. 目标

当前 `computeMonoStatic_HSB` 已经可以先计算有限个 HSB 准确采样点，再对用户指定的单条角度扫描线插值得到 `_RCS_HSB.txt`。如果要同时支持单线扫描和全角域扫描，建议不要再增加额外的“全角域开关”，而是把建库后的输出阶段统一改成由 `scaThetaStart/scaThetaEnd/scaPhiStart/scaPhiEnd/scaStep` 控制的角域网格输出。

目标输出为：

1. 用户配置角域网格上的总 RCS、V 极化 RCS、H 极化 RCS。
2. 用户配置角域网格上的投影复极化电场 `E_v` 和 `E_h`。
3. 输出 HSB 准确采样点信息，用于检查建库点分布、准确点投影复场和插值误差。

## 2. 当前代码状态

HSB 准确点建库位置：

```text
Code/RCS.h
```

函数入口：

```cpp
template<typename SolverType, typename ComputeVFunc>
inline void computeMonoStatic_HSB(SolverType& solver, const RCSExportConfig& cfg,
    ComputeVFunc computeV, double rho = 1.2)
```

当前流程为：

1. 根据包围球半径和 `rho` 计算空间带宽 `W`。
2. 建立 HSB 二维准确采样点库 `rings`。
3. 在每个准确采样点调用 `computeV` 和 `solver.matrix_solver`。
4. 存储该方向的复远场向量 `RCS_pre` 到 `sample.field`。
5. 对用户配置的角域网格进行插值输出。

关键点是：`sample.field` 已经保存了复远场三分量，所以无论输出单线扫描还是全角域扫描，都不需要重新求解，只需要在目标输出角度上调用 `my_hsb::interpolateHSBField(...)`。

## 3. 推荐接口修改

在 `RCSExportConfig` 中只增加输出内容控制项：

```cpp
bool exportHSBComplexField{ false };
bool exportHSBSamples{ false };
```

含义：

- `exportHSBComplexField`: 是否额外输出投影复极化电场 `E_v/E_h`。
- `exportHSBSamples`: 是否输出 HSB 准确采样点信息。

不再新增全角域开关。是否是单线扫描、单点扫描还是全角域扫描，完全由已有单站扫描配置决定：

- `cfg.scaThetaStart`
- `cfg.scaThetaEnd`
- `cfg.scaPhiStart`
- `cfg.scaPhiEnd`
- `cfg.scaStep`

如果要输出完整球面角域，配置为：

```text
scaThetaStart = 0
scaThetaEnd   = 180
scaPhiStart   = 0
scaPhiEnd     = 360
scaStep       = 目标角度步长
```

其中 `phi = 360` 必须输出，用于绘图时闭合曲线或闭合曲面边界。

典型配置含义：

```text
theta 固定、phi 变化：方位向单线扫描
phi 固定、theta 变化：俯仰向单线扫描
theta 变化、phi 变化：二维角域网格扫描
theta 固定、phi 固定：单点输出
```

## 4. 输出文件设计

建议将 RCS 和复电场拆成两个文件，避免一个文件列数过多。

### 4.1 RCS 文件

文件名：

```text
*_RCS_HSB_*.txt
```

列格式：

```text
theta_deg    phi_deg    theta_idx    phi_idx    rcs_db    rcs_v_db    rcs_h_db
```

其中：

```cpp
rcs   = 10.0 * log10((norm(E_v) + norm(E_h)) / (4.0 * Pi));
rcs_v = 10.0 * log10(norm(E_v) / (4.0 * Pi));
rcs_h = 10.0 * log10(norm(E_h) / (4.0 * Pi));
```

### 4.2 复电场文件

文件名：

```text
*_FIELD_HSB_*.txt
```

建议列格式：

```text
theta_deg phi_deg theta_idx phi_idx
Ev_re Ev_im Eh_re Eh_im
```

其中：

- `Ev/Eh` 是用当前观察方向的 `vTh/vPh` 对插值后的 `RCS_pre[0..2]` 投影得到的复极化电场。

### 4.3 HSB 准确采样点文件

文件名：

```text
*_FIELD_HSB_SAMPLES_*.txt
```

建议列格式：

```text
sample_idx ring_idx ring_sample_idx theta_deg phi_deg phi_order
Ev_re Ev_im Eh_re Eh_im
rcs_db rcs_v_db rcs_h_db
```

其中：

- `sample_idx` 是全局准确点序号，从 0 或 1 开始均可，建议文件里保持一致。
- `ring_idx` 对应第几个 theta ring。
- `ring_sample_idx` 对应该 ring 内第几个 phi 点。
- `phi_order` 是该 theta ring 使用的方位向带宽阶数。
- `Ev/Eh` 和 RCS 是在该准确点方向上由准确求解得到的 `RCS_pre[0..2]` 投影得到。
- `RCS_pre[0..2]` 只作为内部计算量，不写入准确点文件。

## 5. 函数结构修改方案

不要再保留“单线输出循环”和“全角域输出循环”两套逻辑。建议在 `RCS.h` 中拆出一个统一的输出辅助函数：

```cpp
template<typename SolverType>
inline void exportHSBScanGrid(
    SolverType& solver,
    const RCSExportConfig& cfg,
    const std::vector<std::vector<HSBMonoSample>>& rings)
```

这个函数根据 `cfg.scaThetaStart/scaThetaEnd/scaPhiStart/scaPhiEnd/scaStep` 自动判断输出网格：

- theta 固定、phi 变化：输出方位向单线。
- phi 固定、theta 变化：输出俯仰向单线。
- theta 和 phi 都变化：输出二维角域网格。
- theta 和 phi 都固定：输出单点。

再增加一个准确点输出辅助函数：

```cpp
template<typename SolverType>
inline void exportHSBSamples(
    SolverType& solver,
    const RCSExportConfig& cfg,
    const std::vector<std::vector<HSBMonoSample>>& rings)
```

`computeMonoStatic_HSB` 建库完成后增加分支：

```cpp
if (cfg.exportHSBSamples) {
    RCSUtils::exportHSBSamples(solver, cfg, rings);
}

RCSUtils::exportHSBScanGrid(solver, cfg, rings);
return;
```

这样有两个好处：

1. 单线扫描和全角域扫描使用同一套插值和输出逻辑，不需要额外开关。
2. 角域输出逻辑和准确点输出逻辑集中在独立函数里，后续修改列格式不会影响准确点建库。

## 6. 统一角域输出循环

统一角域输出循环直接使用单站扫描设置。它不需要知道用户想做“单线”还是“全角域”，只需要按 theta 和 phi 的起止范围生成输出点。

theta 覆盖：

```text
theta = scaThetaStart, scaThetaStart + scaStep, ..., scaThetaEnd
```

phi 覆盖：

```text
phi = scaPhiStart, scaPhiStart + scaStep, ..., scaPhiEnd
```

当 `scaPhiEnd = 360` 时，`phi = 360` 也要输出，不做去重。虽然它和 `phi = 0` 是同一物理方位，但保留该点更方便后处理绘图。

伪代码：

```cpp
auto rcsFile = cfg.openFile("_RCS_HSB", "txt");
std::ofstream fieldFile;
if (cfg.exportHSBComplexField) {
    fieldFile = cfg.openFile("_FIELD_HSB", "txt");
}

auto angleCount = [](double start, double end, double step) {
    if (std::abs(end - start) < 1.0e-9) {
        return 1;
    }
    return static_cast<int>(
        std::floor(std::abs(end - start) / step + 1.0e-9)) + 1;
};

auto angleAt = [](double start, double end, double step, int idx, int count) {
    if (idx == count - 1) {
        return end;
    }
    const double direction = (end >= start) ? 1.0 : -1.0;
    return start + direction * idx * step;
};

const int thetaNum = angleCount(
    cfg.scaThetaStart, cfg.scaThetaEnd, cfg.scaStep);
const int phiNum = angleCount(
    cfg.scaPhiStart, cfg.scaPhiEnd, cfg.scaStep);

for (int it = 0; it < thetaNum; ++it) {
    double thetaDeg = angleAt(
        cfg.scaThetaStart, cfg.scaThetaEnd, cfg.scaStep, it, thetaNum);

    for (int ip = 0; ip < phiNum; ++ip) {
        double phiDeg = angleAt(
            cfg.scaPhiStart, cfg.scaPhiEnd, cfg.scaStep, ip, phiNum);

        const double phiInterpDeg = my_hsb::normalizePhiDeg(phiDeg);

        std::vector<std::complex<double>> RCS_pre =
            my_hsb::interpolateHSBField(rings, thetaDeg, phiInterpDeg);

        EMSource pw;
        pw.initPW(data_pol, thetaDeg, phiDeg);

        double vec_v[3] = { pw.vTh(0), pw.vTh(1), pw.vTh(2) };
        double vec_h[3] = { pw.vPh(0), pw.vPh(1), pw.vPh(2) };

        std::complex<double> E_v =
            vec_v[0] * RCS_pre[0] +
            vec_v[1] * RCS_pre[1] +
            vec_v[2] * RCS_pre[2];

        std::complex<double> E_h =
            vec_h[0] * RCS_pre[0] +
            vec_h[1] * RCS_pre[1] +
            vec_h[2] * RCS_pre[2];

        double rcs = 10.0 * std::log10(
            (std::norm(E_v) + std::norm(E_h)) / (4.0 * Pi));
        double rcs_v = 10.0 * std::log10(std::norm(E_v) / (4.0 * Pi));
        double rcs_h = 10.0 * std::log10(std::norm(E_h) / (4.0 * Pi));

        rcsFile << thetaDeg << "\t" << phiDeg << "\t"
            << it << "\t" << ip << "\t"
            << rcs << "\t" << rcs_v << "\t" << rcs_h << "\n";

        if (cfg.exportHSBComplexField) {
            fieldFile << thetaDeg << "\t" << phiDeg << "\t"
                << it << "\t" << ip << "\t"
                << E_v.real() << "\t" << E_v.imag() << "\t"
                << E_h.real() << "\t" << E_h.imag() << "\n";
        }
    }
}
```

注意：

1. 上面的伪代码中需要把 `data_pol` 传入辅助函数，或者在辅助函数内部重新通过 `cfg.polWaveStr()` 计算。
2. 文件中写出的 `phiDeg` 保留原始配置角度，所以当 `scaPhiEnd = 360` 时文件会写出 `360`。
3. 插值时使用 `phiInterpDeg = normalizePhiDeg(phiDeg)`，避免 HSB 插值内部因为 `360` 和 `0` 的周期等价关系出错。

## 7. HSB 准确点信息输出方案

角域输出是插值场。为了检查 HSB 建库是否合理、插值误差是否可控，建议将准确采样点输出作为独立功能，而不是只在调试时临时打印。

### 7.1 输出时机

准确点信息应在 `rings` 完成后统一输出，而不是在每个准确点求解时直接写文件。

推荐位置：

```cpp
// HSB accurate samples have been solved and stored in rings.
if (cfg.exportHSBSamples) {
    RCSUtils::exportHSBSamples(solver, cfg, rings);
}
```

原因：

1. `rings` 已经包含所有准确点的角度、`phiOrder` 和内部复远场，可用于计算 `E_v/E_h` 和 RCS。
2. 输出逻辑不影响 `matrix_solver` 主循环。
3. 后续可以复用同一份 `rings` 同时输出准确点和用户配置角域的插值结果。

### 7.2 准确点输出伪代码

```cpp
template<typename SolverType>
inline void exportHSBSamples(
    SolverType& solver,
    const RCSExportConfig& cfg,
    const std::vector<std::vector<HSBMonoSample>>& rings)
{
    auto sampleFile = cfg.openFile("_FIELD_HSB_SAMPLES", "txt");
    const std::string pol_wave = cfg.polWaveStr();
    const double data_pol = (pol_wave == "V") ? 0.0 : 90.0;

    sampleFile << "sample_idx\tring_idx\tring_sample_idx\t"
        << "theta_deg\tphi_deg\tphi_order\t"
        << "Ev_re\tEv_im\tEh_re\tEh_im\t"
        << "rcs_db\trcs_v_db\trcs_h_db\n";

    int sampleIdx = 0;
    for (int ringIdx = 0; ringIdx < static_cast<int>(rings.size()); ++ringIdx) {
        const auto& ring = rings[ringIdx];

        for (int ringSampleIdx = 0;
             ringSampleIdx < static_cast<int>(ring.size());
             ++ringSampleIdx) {
            const HSBMonoSample& sample = ring[ringSampleIdx];
            const std::vector<std::complex<double>>& RCS_pre = sample.field;

            EMSource pw;
            pw.initPW(data_pol, sample.thetaDeg, sample.phiDeg);

            double vec_v[3] = { pw.vTh(0), pw.vTh(1), pw.vTh(2) };
            double vec_h[3] = { pw.vPh(0), pw.vPh(1), pw.vPh(2) };

            std::complex<double> E_v =
                vec_v[0] * RCS_pre[0] +
                vec_v[1] * RCS_pre[1] +
                vec_v[2] * RCS_pre[2];

            std::complex<double> E_h =
                vec_h[0] * RCS_pre[0] +
                vec_h[1] * RCS_pre[1] +
                vec_h[2] * RCS_pre[2];

            double rcs = 10.0 * std::log10(
                (std::norm(E_v) + std::norm(E_h)) / (4.0 * Pi));
            double rcs_v = 10.0 * std::log10(std::norm(E_v) / (4.0 * Pi));
            double rcs_h = 10.0 * std::log10(std::norm(E_h) / (4.0 * Pi));

            sampleFile << sampleIdx << "\t"
                << ringIdx << "\t" << ringSampleIdx << "\t"
                << sample.thetaDeg << "\t" << sample.phiDeg << "\t"
                << sample.phiOrder << "\t"
                << E_v.real() << "\t" << E_v.imag() << "\t"
                << E_h.real() << "\t" << E_h.imag() << "\t"
                << rcs << "\t" << rcs_v << "\t" << rcs_h << "\n";

            ++sampleIdx;
        }
    }
}
```

### 7.3 准确点输出和进度输出的区别

准确点文件输出和控制台进度输出是两个不同目的：

- 控制台输出用于运行时观察当前算到第几个点。
- `_FIELD_HSB_SAMPLES` 文件用于保存全部准确点的数值结果，便于后处理和误差分析。

建议两者都保留。控制台输出放在 `solver.matrix_solver(...)` 后；准确点文件输出放在建库完成后。

## 8. 扫描形态由角度设置控制

不再使用额外开关区分单线扫描和全角域扫描。`computeMonoStatic_HSB` 建库完成后总是调用统一输出函数：

```cpp
RCSUtils::exportHSBScanGrid(solver, cfg, rings);
```

扫描形态由配置自然决定：

```text
theta 固定、phi 变化：方位向单线扫描
phi 固定、theta 变化：俯仰向单线扫描
theta 变化、phi 变化：二维角域扫描
theta 固定、phi 固定：单点输出
```

准确点输出不应该独占流程。也就是说：

```cpp
if (cfg.exportHSBSamples) {
    exportHSBSamples(...);
}
```

执行后继续调用 `exportHSBScanGrid(...)` 输出用户配置的角域。

## 9. 推荐实现顺序

1. 在 `RCSExportConfig.h` 增加复场输出和准确点输出配置项。
2. 在 `RCS.h` 增加 `exportHSBSamples(...)` 辅助函数。
3. 在 `computeMonoStatic_HSB` 准确点建库完成后添加 `cfg.exportHSBSamples` 分支。
4. 先输出 `_FIELD_HSB_SAMPLES.txt` 并编译，确认准确点数量与控制台 `HSB sample points` 一致。
5. 在 `RCS.h` 增加 `exportHSBScanGrid(...)` 辅助函数。
6. 用 `exportHSBScanGrid(...)` 替换 `computeMonoStatic_HSB` 末尾原有的单线输出循环。
7. 先只输出 `_RCS_HSB.txt` 并编译。
8. 再增加 `_FIELD_HSB.txt` 复场输出。
9. 测试单线配置，例如 `scaThetaStart=90`、`scaThetaEnd=90`、`scaPhiStart=0`、`scaPhiEnd=360`、`scaStep=10`。
10. 测试二维配置，例如 `scaThetaStart=0`、`scaThetaEnd=180`、`scaPhiStart=0`、`scaPhiEnd=360`、`scaStep=10`。
11. 再切换到目标精度，例如 `scaStep=1` 或更细。

## 10. 验证建议

### 10.1 点数验证

例如：

```text
scaThetaStart = 0
scaThetaEnd   = 180
scaPhiStart   = 0
scaPhiEnd     = 360
scaStep       = 1
```

则输出点数应为：

```text
181 * 361 = 65341
```

这里包含 `phi = 360`，用于绘图闭合。

### 10.2 数值验证

选取几条典型切线：

- `theta = 90, phi = 0..360`
- `phi = 0, theta = 0..180`
- `phi = 90, theta = 0..180`

分别与原始逐角度 `computeMonoStatic` 或 `computeMonoStatic_PMCHWT` 结果对比。

### 10.3 文件验证

检查：

- RCS 文件是否每行都有 `theta/phi/rcs`。
- 复场文件是否每行都有 `Ev/Eh` 的实部和虚部。
- 准确点文件行数是否等于 `HSB sample points`。
- 准确点文件中的 `ring_idx/ring_sample_idx/phi_order` 是否能还原 HSB 采样结构。
- 当 `scaPhiEnd = 360` 时，输出文件是否包含 `phi = 360`。
- `theta = 0` 和 `theta = 180` 极点附近的 `E_v/E_h` 是否稳定。

## 11. 注意事项

1. 不建议插值 dB RCS。当前代码插值的是复远场 `RCS_pre`，这是正确方向。
2. 全角域输出点数可能很大，但它只做插值和写文件，不会重新调用 `matrix_solver`。
3. 如果输出 `1` 度全角域，复场文件会有六万多行，文件体积会明显增加。
4. 极点 `theta = 0/180` 的水平/垂直基向量定义可能存在方向约定问题。RCS 总量通常更稳定，极化分量需要结合 `EMSource::initPW` 的定义解释。
5. 如果后续要做极化方向图，建议优先使用 `E_v/E_h` 的复数结果，而不是只用 dB RCS。
