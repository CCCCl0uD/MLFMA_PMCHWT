# lexpDIE 抽取修改方案

## 目标

将 `Code/MLFMM.cpp` 中 PMCHWT/DIE 分支里的 finest-level disaggregation 和 near-field 汇聚代码抽取为 `ProcessMLFMM::lexpDIE`，放到 `Code/ProcessMLFMM.h` 中，形成和 `lexpPEC` 对应的介质体/介质界面方程后处理函数。

抽取后的 `MLFMM::applyAx` PMCHWT 分支应只保留：

1. 四次远场传播：
   - `J1`: `mexp -> m2m -> m2l -> l2l`
   - `M1`: `mexp -> m2m -> m2l -> l2l`
   - `J2`: `mexp -> m2m -> m2l -> l2l`
   - `M2`: `mexp -> m2m -> m2l -> l2l`
2. 一次 `ProcessMLFMM::lexpDIE(...)`，负责把四组 `Bm` 转换为矩阵向量乘结果 `w`。

## 当前代码结构

当前 PMCHWT 分支中，远场 disaggregation 和近场汇聚直接写在 `MLFMM::applyAx` 内：

```cpp
auto integrateL = [](const std::vector<Complex3D>& V,
    const std::vector<Complex3D>& B) -> std::complex<double> {
    ...
};

auto integrateK = [](const std::vector<Complex3D>& V,
    const std::vector<Complex3D>& B,
    const std::vector<kp_Point>& kp) -> std::complex<double> {
    ...
};

#pragma omp parallel for schedule(dynamic)
for (int nodeIdx = 0; nodeIdx < finestN; ++nodeIdx) {
    ...
    w[rwgid] += const1 * L1J;
    w[rwgid] += -const1 * K1M;
    w[N + rwgid] += const1 * K1J;
    w[N + rwgid] += const1 * L1M;
    w[rwgid] += const2 * etai * L2J;
    w[rwgid] += -const2 * K2M;
    w[N + rwgid] += const2 * K2J;
    w[N + rwgid] += const2 * L2M / etai;
}

#pragma omp parallel for schedule(dynamic)
for (int nodeIdx = 0; nodeIdx < finestN; ++nodeIdx) {
    ...
    w[rwgid] += near-field block;
    w[N + rwgid] += near-field block;
}
```

这段代码逻辑完整，但放在 `applyAx` 里会让 PMCHWT 分支过长，也和 `ProcessMLFMM::lexpPEC` 的职责不一致。

## 建议的函数位置

在 `Code/ProcessMLFMM.h` 中，放在 `lexpPEC` 后面，新增：

```cpp
inline void lexpDIE(...);
```

建议顺序：

1. `lexpKpDot`
2. `lexpPEC`
3. `lexpDIE`

如果后续希望 `lexpPEC` 和 `lexpDIE` 共用小工具，可在 `lexpKpDot` 附近再新增：

```cpp
inline std::complex<double> lexpLDot(...);
inline std::complex<double> lexpKDot(...);
```

本次可以先把 `integrateL/integrateK` 作为 `lexpDIE` 内部 lambda，改动更集中。

## `lexpDIE` 建议接口

推荐接口如下：

```cpp
inline void lexpDIE(std::complex<double>* w,
    const std::vector<OCTree::Node*>& finestNodes,
    const std::vector<std::vector<Complex3D>>& Vfmj1,
    const std::vector<std::vector<Complex3D>>& Vfmj2,
    const std::vector<std::vector<std::vector<Complex3D>>>& Bm_J1,
    const std::vector<std::vector<std::vector<Complex3D>>>& Bm_M1,
    const std::vector<std::vector<std::vector<Complex3D>>>& Bm_J2,
    const std::vector<std::vector<std::vector<Complex3D>>>& Bm_M2,
    const std::vector<kp_Point>& kp_k1,
    const std::vector<kp_Point>& kp_k2,
    const std::vector<double>& WGL_k1,
    double WGL_phi_k1,
    const std::vector<double>& WGL_k2,
    double WGL_phi_k2,
    int maxLevel_,
    int N,
    std::complex<double> const1,
    std::complex<double> const2,
    std::complex<double> etai,
    const std::vector<std::vector<std::complex<double>>>& Z_near,
    const std::vector<std::vector<int>>& Z_near_id,
    const std::complex<double>* x)
```

参数含义：

- `w`: 输出向量，长度应为 `2 * N`。
- `finestNodes`: `octreeNodes_[maxlvl]`。
- `Vfmj1`: k1 对应的测试函数角谱，即当前 `Vfmj`。
- `Vfmj2`: k2 对应的测试函数角谱，即当前 `Vfmj2`。
- `Bm_J1/Bm_M1/Bm_J2/Bm_M2`: 四次远场传播得到的局部展开。
- `kp_k1/kp_k2`: 最细层 k1/k2 方向采样点，即 `kp_lvl_k1[0]` 和 `kp_lvl_k2[0]`。
- `WGL_k1/WGL_phi_k1`: 最细层 k1 积分权重。
- `WGL_k2/WGL_phi_k2`: 最细层 k2 积分权重。
- `maxLevel_`: 用于定位 finest-level Bm 索引 `maxLevel_ - 2`。
- `N`: 单个未知量块大小，即 `row`。
- `const1`: `wave.k1() * wave.k1() / (16.0 * Pi * Pi)`。
- `const2`: `wave.k2() * wave.k2() / (16.0 * Pi * Pi)`。
- `etai`: `wave.etai()`。
- `Z_near/Z_near_id/x`: 近场 CSR 数据和输入向量。

## 权重处理建议

由于当前 `m2l` 已经不再乘积分权重，`lexpDIE` 应和 `lexpPEC` 保持一致，在最终角谱积分累加时乘入权重。

因此建议把原来的：

```cpp
val += dot(V[idx], B[idx]);
```

改为：

```cpp
const int t = idx / pN;
const double weight = WGL[t] * WGL_phi;
val += dot(V[idx], B[idx]) * weight;
```

`integrateK` 同理：

```cpp
val += dot(kCrossV, B[idx]) * weight;
```

为此，`lexpDIE` 内部需要知道 `pN`，可通过：

```cpp
const int pN_k1 = static_cast<int>(kp_k1.size() / WGL_k1.size());
const int pN_k2 = static_cast<int>(kp_k2.size() / WGL_k2.size());
```

得到。为了更稳妥，也可以把 `pN_k1`、`pN_k2` 显式作为参数传入。推荐显式传入，避免在权重尺寸异常时静默出错：

```cpp
int pN_k1,
int pN_k2,
```

最终接口可在 `kp_k1/kp_k2` 后加入这两个参数。

## `lexpDIE` 内部实现结构

建议函数内部分成三段：

### 1. 局部积分 lambda

```cpp
auto integrateL = [](const std::vector<Complex3D>& V,
    const std::vector<Complex3D>& B,
    const std::vector<double>& WGL,
    double WGL_phi,
    int pN) -> std::complex<double> {
    std::complex<double> val(0.0, 0.0);
    const int K = static_cast<int>(V.size());
    for (int idx = 0; idx < K; ++idx) {
        const int t = idx / pN;
        const double weight = WGL[t] * WGL_phi;
        val += dot(V[idx], B[idx]) * weight;
    }
    return val;
};
```

```cpp
auto integrateK = [](const std::vector<Complex3D>& V,
    const std::vector<Complex3D>& B,
    const std::vector<kp_Point>& kp,
    const std::vector<double>& WGL,
    double WGL_phi,
    int pN) -> std::complex<double> {
    std::complex<double> val(0.0, 0.0);
    const int K = static_cast<int>(V.size());
    for (int idx = 0; idx < K; ++idx) {
        const double khat[3] = { -kp[idx].x, -kp[idx].y, -kp[idx].z };

        Complex3D kCrossV;
        cross(kCrossV, khat, V[idx]);

        const int t = idx / pN;
        const double weight = WGL[t] * WGL_phi;
        val += dot(kCrossV, B[idx]) * weight;
    }
    return val;
};
```

### 2. 远场汇聚

把当前 `applyAx` 中的 finest-level disaggregation 循环完整搬入函数：

```cpp
const int finestIdx = maxLevel_ - 2;
const int finestN = static_cast<int>(finestNodes.size());

#pragma omp parallel for schedule(dynamic)
for (int nodeIdx = 0; nodeIdx < finestN; ++nodeIdx) {
    OCTree::Node* nd = finestNodes[nodeIdx];

    const auto& BJ1 = Bm_J1[finestIdx][nodeIdx];
    const auto& BM1 = Bm_M1[finestIdx][nodeIdx];
    const auto& BJ2 = Bm_J2[finestIdx][nodeIdx];
    const auto& BM2 = Bm_M2[finestIdx][nodeIdx];

    for (int ri = 0; ri < static_cast<int>(nd->rwgIndices.size()); ++ri) {
        const int rwgid = nd->rwgIndices[ri]->rwgid;

        const auto& V1 = Vfmj1[rwgid];
        const auto& V2 = Vfmj2[rwgid];

        const auto L1J = integrateL(V1, BJ1, WGL_k1, WGL_phi_k1, pN_k1);
        const auto K1J = integrateK(V1, BJ1, kp_k1, WGL_k1, WGL_phi_k1, pN_k1);
        const auto K1M = integrateK(V1, BM1, kp_k1, WGL_k1, WGL_phi_k1, pN_k1);
        const auto L1M = integrateL(V1, BM1, WGL_k1, WGL_phi_k1, pN_k1);

        const auto L2J = integrateL(V2, BJ2, WGL_k2, WGL_phi_k2, pN_k2);
        const auto K2J = integrateK(V2, BJ2, kp_k2, WGL_k2, WGL_phi_k2, pN_k2);
        const auto K2M = integrateK(V2, BM2, kp_k2, WGL_k2, WGL_phi_k2, pN_k2);
        const auto L2M = integrateL(V2, BM2, WGL_k2, WGL_phi_k2, pN_k2);

        w[rwgid] += const1 * L1J;
        w[rwgid] += -const1 * K1M;

        w[N + rwgid] += const1 * K1J;
        w[N + rwgid] += const1 * L1M;

        w[rwgid] += const2 * etai * L2J;
        w[rwgid] += -const2 * K2M;

        w[N + rwgid] += const2 * K2J;
        w[N + rwgid] += const2 * L2M / etai;
    }
}
```

### 3. 近场汇聚

把当前 near-field 循环也搬入 `lexpDIE`：

```cpp
#pragma omp parallel for schedule(dynamic)
for (int nodeIdx = 0; nodeIdx < finestN; ++nodeIdx) {
    OCTree::Node* nd = finestNodes[nodeIdx];
    for (int ri = 0; ri < static_cast<int>(nd->rwgIndices.size()); ++ri) {
        const int rwgid = nd->rwgIndices[ri]->rwgid;

        for (int j = 0; j < static_cast<int>(Z_near[rwgid].size()); ++j) {
            w[rwgid] += x[Z_near_id[rwgid][j]] * Z_near[rwgid][j];
        }

        for (int j = 0; j < static_cast<int>(Z_near[N + rwgid].size()); ++j) {
            w[N + rwgid] += x[Z_near_id[N + rwgid][j]] * Z_near[N + rwgid][j];
        }
    }
}
```

## `MLFMM::applyAx` 调用方式

PMCHWT 分支中，四次 `l2l` 后，用一行调用替换当前的远场和近场汇聚代码：

```cpp
const int pN_k1 = static_cast<int>(phi_level_k1[0].size());
const int pN_k2 = static_cast<int>(phi_level_k2[0].size());

ProcessMLFMM::lexpDIE(w,
    octreeNodes_[maxlvl],
    Vfmj,
    Vfmj2,
    Bm_J1,
    Bm_M1,
    Bm_J2,
    Bm_M2,
    kp_lvl_k1[0],
    kp_lvl_k2[0],
    WGL_k1[0],
    WGL_phi_k1[0],
    WGL_k2[0],
    WGL_phi_k2[0],
    pN_k1,
    pN_k2,
    maxlvl,
    N,
    const1,
    const2,
    etai,
    Z_near,
    Z_near_id,
    x);
```

然后删除 `applyAx` 内原来的：

- `integrateL` lambda
- `integrateK` lambda
- far-field finest-level disaggregation 循环
- near-field 循环

## 建议的最终接口版本

综合权重和 `pN` 的显式传入，建议最终采用：

```cpp
inline void lexpDIE(std::complex<double>* w,
    const std::vector<OCTree::Node*>& finestNodes,
    const std::vector<std::vector<Complex3D>>& Vfmj1,
    const std::vector<std::vector<Complex3D>>& Vfmj2,
    const std::vector<std::vector<std::vector<Complex3D>>>& Bm_J1,
    const std::vector<std::vector<std::vector<Complex3D>>>& Bm_M1,
    const std::vector<std::vector<std::vector<Complex3D>>>& Bm_J2,
    const std::vector<std::vector<std::vector<Complex3D>>>& Bm_M2,
    const std::vector<kp_Point>& kp_k1,
    const std::vector<kp_Point>& kp_k2,
    const std::vector<double>& WGL_k1,
    double WGL_phi_k1,
    const std::vector<double>& WGL_k2,
    double WGL_phi_k2,
    int pN_k1,
    int pN_k2,
    int maxLevel_,
    int N,
    std::complex<double> const1,
    std::complex<double> const2,
    std::complex<double> etai,
    const std::vector<std::vector<std::complex<double>>>& Z_near,
    const std::vector<std::vector<int>>& Z_near_id,
    const std::complex<double>* x)
```

## 并行安全性

当前两段循环都以 `nodeIdx` 并行。每个 RWG 应只属于一个 finest node，因此每个线程写入的 `w[rwgid]` 和 `w[N + rwgid]` 不应重复。

抽取为 `lexpDIE` 后仍保持相同的并行结构，不改变并行写入模式。

需要注意：如果将来一个 RWG 可能出现在多个 finest node 中，就需要对 `w` 写入改为归约、原子操作或按 node 暂存后合并。本次不改变这个假设。

## 验证方案

1. 编译 `MLFMA.sln`，确认 `ProcessMLFMM::lexpDIE` 的声明、定义和调用参数完全匹配。
2. 运行 PMCHWT 模式，确认 `applyAx` 可正常完成。
3. 对比抽取前后的 `matvec_mlfma.txt` 或同一输入下的残差历史，确认数值变化只来自预期的权重位置调整。
4. 全仓搜索确认 `MLFMM::applyAx` 中不再保留 PMCHWT 专用的 `integrateL/integrateK` lambda。
5. 全仓搜索确认 DIE/PMCHWT 的远场最终积分权重只在 `lexpDIE` 中乘入，而不是回到 `m2l`。

## 预期收益

- `MLFMM::applyAx` 的 PMCHWT 分支会更短，只负责调度多层传播流程。
- `ProcessMLFMM.h` 会集中管理 PEC 和 DIE 两类 finest-level disaggregation。
- 权重位置和 `lexpPEC` 保持一致，后续检查“是否重复乘权重”会更容易。
- 后续如果要优化 `integrateL/integrateK`，只需要改 `lexpDIE` 内部，不必进入 `MLFMM::applyAx` 主流程。

## 完整代码草案

下面给出可以直接写入代码的完整版本。代码注释刻意使用 ASCII，避免再次出现中文编码乱码。

### 1. 写入 `Code/ProcessMLFMM.h` 的完整 `lexpDIE`

建议放在 `lexpPEC` 后面、`namespace ProcessMLFMM` 结束之前。

```cpp
// ============================================================================
// 8. lexpDIE - DIE/PMCHWT finest-level disaggregation + near field
// ============================================================================
inline void lexpDIE(std::complex<double>* w,
    const std::vector<OCTree::Node*>& finestNodes,
    const std::vector<std::vector<Complex3D>>& Vfmj1,
    const std::vector<std::vector<Complex3D>>& Vfmj2,
    const std::vector<std::vector<std::vector<Complex3D>>>& Bm_J1,
    const std::vector<std::vector<std::vector<Complex3D>>>& Bm_M1,
    const std::vector<std::vector<std::vector<Complex3D>>>& Bm_J2,
    const std::vector<std::vector<std::vector<Complex3D>>>& Bm_M2,
    const std::vector<kp_Point>& kp_k1,
    const std::vector<kp_Point>& kp_k2,
    const std::vector<double>& WGL_k1,
    double WGL_phi_k1,
    const std::vector<double>& WGL_k2,
    double WGL_phi_k2,
    int pN_k1,
    int pN_k2,
    int maxLevel_,
    int N,
    std::complex<double> const1,
    std::complex<double> const2,
    std::complex<double> etai,
    const std::vector<std::vector<std::complex<double>>>& Z_near,
    const std::vector<std::vector<int>>& Z_near_id,
    const std::complex<double>* x)
{
    auto integrateL = [](const std::vector<Complex3D>& V,
        const std::vector<Complex3D>& B,
        const std::vector<double>& WGL,
        double WGL_phi,
        int pN) -> std::complex<double> {
        std::complex<double> val(0.0, 0.0);
        const int K = static_cast<int>(V.size());
        for (int idx = 0; idx < K; ++idx) {
            const int t = idx / pN;
            const double weight = WGL[t] * WGL_phi;
            val += dot(V[idx], B[idx]) * weight;
        }
        return val;
    };

    auto integrateK = [](const std::vector<Complex3D>& V,
        const std::vector<Complex3D>& B,
        const std::vector<kp_Point>& kp,
        const std::vector<double>& WGL,
        double WGL_phi,
        int pN) -> std::complex<double> {
        std::complex<double> val(0.0, 0.0);
        const int K = static_cast<int>(V.size());
        for (int idx = 0; idx < K; ++idx) {
            const double khat[3] = { -kp[idx].x, -kp[idx].y, -kp[idx].z };

            Complex3D kCrossV;
            cross(kCrossV, khat, V[idx]);

            const int t = idx / pN;
            const double weight = WGL[t] * WGL_phi;
            val += dot(kCrossV, B[idx]) * weight;
        }
        return val;
    };

    const int finestIdx = maxLevel_ - 2;
    const int finestN = static_cast<int>(finestNodes.size());

#pragma omp parallel for schedule(dynamic)
    for (int nodeIdx = 0; nodeIdx < finestN; ++nodeIdx) {
        OCTree::Node* nd = finestNodes[nodeIdx];

        const auto& BJ1 = Bm_J1[finestIdx][nodeIdx];
        const auto& BM1 = Bm_M1[finestIdx][nodeIdx];
        const auto& BJ2 = Bm_J2[finestIdx][nodeIdx];
        const auto& BM2 = Bm_M2[finestIdx][nodeIdx];

        for (int ri = 0; ri < static_cast<int>(nd->rwgIndices.size()); ++ri) {
            const int rwgid = nd->rwgIndices[ri]->rwgid;

            const auto& V1 = Vfmj1[rwgid];
            const auto& V2 = Vfmj2[rwgid];

            const auto L1J = integrateL(V1, BJ1, WGL_k1, WGL_phi_k1, pN_k1);
            const auto K1J = integrateK(V1, BJ1, kp_k1, WGL_k1, WGL_phi_k1, pN_k1);
            const auto K1M = integrateK(V1, BM1, kp_k1, WGL_k1, WGL_phi_k1, pN_k1);
            const auto L1M = integrateL(V1, BM1, WGL_k1, WGL_phi_k1, pN_k1);

            const auto L2J = integrateL(V2, BJ2, WGL_k2, WGL_phi_k2, pN_k2);
            const auto K2J = integrateK(V2, BJ2, kp_k2, WGL_k2, WGL_phi_k2, pN_k2);
            const auto K2M = integrateK(V2, BM2, kp_k2, WGL_k2, WGL_phi_k2, pN_k2);
            const auto L2M = integrateL(V2, BM2, WGL_k2, WGL_phi_k2, pN_k2);

            w[rwgid] += const1 * L1J;
            w[rwgid] += -const1 * K1M;

            w[N + rwgid] += const1 * K1J;
            w[N + rwgid] += const1 * L1M;

            w[rwgid] += const2 * etai * L2J;
            w[rwgid] += -const2 * K2M;

            w[N + rwgid] += const2 * K2J;
            w[N + rwgid] += const2 * L2M / etai;
        }
    }

#pragma omp parallel for schedule(dynamic)
    for (int nodeIdx = 0; nodeIdx < finestN; ++nodeIdx) {
        OCTree::Node* nd = finestNodes[nodeIdx];
        for (int ri = 0; ri < static_cast<int>(nd->rwgIndices.size()); ++ri) {
            const int rwgid = nd->rwgIndices[ri]->rwgid;

            for (int j = 0; j < static_cast<int>(Z_near[rwgid].size()); ++j) {
                w[rwgid] += x[Z_near_id[rwgid][j]] * Z_near[rwgid][j];
            }

            for (int j = 0; j < static_cast<int>(Z_near[N + rwgid].size()); ++j) {
                w[N + rwgid] += x[Z_near_id[N + rwgid][j]] * Z_near[N + rwgid][j];
            }
        }
    }
}
```

### 2. `Code/MLFMM.cpp` 中 PMCHWT 分支的替换代码

四次 `mexp -> m2m -> m2l -> l2l` 保留不变。把原来的 `integrateL`、`integrateK`、far-field 循环和 near-field 循环整体删除，替换为：

```cpp
const int pN_k1 = static_cast<int>(phi_level_k1[0].size());
const int pN_k2 = static_cast<int>(phi_level_k2[0].size());

ProcessMLFMM::lexpDIE(w,
    octreeNodes_[maxlvl],
    Vfmj,
    Vfmj2,
    Bm_J1,
    Bm_M1,
    Bm_J2,
    Bm_M2,
    kp_lvl_k1[0],
    kp_lvl_k2[0],
    WGL_k1[0],
    WGL_phi_k1[0],
    WGL_k2[0],
    WGL_phi_k2[0],
    pN_k1,
    pN_k2,
    maxlvl,
    N,
    const1,
    const2,
    etai,
    Z_near,
    Z_near_id,
    x);

return;
```

### 3. 替换后的 PMCHWT 分支末尾应类似这样

下面这段展示从四次 `l2l` 结束到 `return` 的完整形态：

```cpp
// ---- Pass 4: k2, M -> K2(M), L2(M) ----
ProcessMLFMM::mexp(Sm_M2, octreeNodes_[maxlvl], Vsmi2, kp_lvl_k2[0], x, N, maxlvl);
ProcessMLFMM::m2m(Sm_M2, octreeNodes_, kp_lvl_k2, interpol_k2, maxlvl, levelSpan, wave.k2());
ProcessMLFMM::m2l(Sm_M2, Bm_M2, octreeNodes_, TF2, kp_lvl_k2, maxlvl, levelSpan);
ProcessMLFMM::l2l(Bm_M2, octreeNodes_, kp_lvl_k2, phi_level_k2, interpol_k2, maxlvl, levelSpan, wave.k2());

const int pN_k1 = static_cast<int>(phi_level_k1[0].size());
const int pN_k2 = static_cast<int>(phi_level_k2[0].size());

ProcessMLFMM::lexpDIE(w,
    octreeNodes_[maxlvl],
    Vfmj,
    Vfmj2,
    Bm_J1,
    Bm_M1,
    Bm_J2,
    Bm_M2,
    kp_lvl_k1[0],
    kp_lvl_k2[0],
    WGL_k1[0],
    WGL_phi_k1[0],
    WGL_k2[0],
    WGL_phi_k2[0],
    pN_k1,
    pN_k2,
    maxlvl,
    N,
    const1,
    const2,
    etai,
    Z_near,
    Z_near_id,
    x);

return;
```

### 4. 可选的防御性检查

如果希望在调试阶段更早发现维度错误，可以在 `lexpDIE` 开头加入：

```cpp
if (pN_k1 <= 0 || pN_k2 <= 0) {
    throw std::invalid_argument("lexpDIE invalid phi count.");
}
if (kp_k1.size() % WGL_k1.size() != 0 || kp_k2.size() % WGL_k2.size() != 0) {
    throw std::invalid_argument("lexpDIE inconsistent kp/WGL size.");
}
```

注意：如果加入 `throw std::invalid_argument`，需要确认当前头文件包含链已经提供 `<stdexcept>`；否则应在合适的公共头文件或 `ProcessMLFMM.h` 中包含它。
