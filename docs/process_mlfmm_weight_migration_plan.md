# ProcessMLFMM 权重迁移修改方案

## 目标

将 `ProcessMLFMM.h` 中 `m2l` 阶段的积分权重移除，使 `m2l` 只负责执行

```cpp
Bm += TF * Sm;
```

积分权重 `WGL[t] * WGL_phi` 只在 `lexpPEC` 的最终积分累加处乘入，也就是只作用于

```cpp
lexpKpDot(B[kp], V[kp])
```

这样可以避免权重在多层传播或不同算子路径中被提前混入 `Bm`。

## 当前状态

`Code/ProcessMLFMM.h` 中 `m2l` 的签名仍接收：

```cpp
const std::vector<std::vector<double>>& WGL,
const std::vector<double>& WGL_phi,
const std::vector<std::vector<double>>& phi_level
```

但函数体内实际已经没有乘权重，原权重代码处于注释状态：

```cpp
// int t = kp_i / pN;
// double weight = WGL[rev][t] * WGL_phi[rev];
// Bmvec[kp_i].x += TFvec[kp_i] * Smvec[kp_i].x * weight;
```

PEC 远场汇聚目前在 `Code/MLFMM.cpp` 的 `MLFMM::applyAx` 内手写循环完成，未调用 `ProcessMLFMM::lexpPEC`。该循环目前也没有乘 `WGL_k1[0][tt] * WGL_phi_k1[0]`。

## 修改方案

### 1. 简化 `m2l` 接口

从 `ProcessMLFMM::m2l` 的形参中删除：

```cpp
const std::vector<std::vector<double>>& WGL,
const std::vector<double>& WGL_phi,
const std::vector<std::vector<double>>& phi_level
```

同时删除函数内部不再需要的：

```cpp
const int pN = static_cast<int>(phi_level[rev].size());
```

以及已经注释掉的旧权重代码。

保留的 `m2l` 核心逻辑应为：

```cpp
Bmvec[kp_i].x += TFvec[kp_i] * Smvec[kp_i].x;
Bmvec[kp_i].y += TFvec[kp_i] * Smvec[kp_i].y;
Bmvec[kp_i].z += TFvec[kp_i] * Smvec[kp_i].z;
```

### 2. 更新所有 `m2l` 调用点

`Code/MLFMM.cpp` 中所有调用需要从：

```cpp
ProcessMLFMM::m2l(Sm, Bm, octreeNodes_, TF,
    kp_lvl_k1, WGL_k1, WGL_phi_k1, phi_level_k1, maxlvl, levelSpan);
```

改为：

```cpp
ProcessMLFMM::m2l(Sm, Bm, octreeNodes_, TF,
    kp_lvl_k1, maxlvl, levelSpan);
```

`k2`、PMCHWT 四个 pass 中的调用同样按此规则修改。

### 3. 扩展 `lexpPEC` 形参

给 `ProcessMLFMM::lexpPEC` 增加最细层积分权重参数：

```cpp
const std::vector<double>& WGL,
double WGL_phi
```

建议放在 `tN, pN` 后面，便于表达它们和最细层角向网格对应：

```cpp
inline void lexpPEC(std::complex<double>* w,
    const std::vector<OCTree::Node*>& finestNodes,
    const std::vector<std::vector<Complex3D>>& Vfmj,
    const std::vector<std::vector<std::vector<Complex3D>>>& Bm,
    int tN, int pN,
    const std::vector<double>& WGL,
    double WGL_phi,
    int maxLevel_,
    std::complex<double> k_, std::complex<double> eta_,
    const std::vector<std::vector<std::complex<double>>>& Z_near,
    const std::vector<std::vector<int>>& Z_near_id,
    const std::complex<double>* x)
```

### 4. 只在 `lexpPEC` 中乘权重

将 `lexpPEC` 内的积分累加从：

```cpp
val += lexpKpDot(B[tt * pN + pp], V[tt * pN + pp]);
```

改为：

```cpp
const int kpIdx = tt * pN + pp;
const double weight = WGL[tt] * WGL_phi;
val += lexpKpDot(B[kpIdx], V[kpIdx]) * weight;
```

这样权重只在最终把角谱积分投影到测试基函数时使用。

### 5. 在 PEC 分支调用 `lexpPEC`

`MLFMM::applyAx` 的 PEC 分支完成：

```cpp
clearSmBm -> mexp -> m2m -> m2l -> l2l
```

后，可用 `ProcessMLFMM::lexpPEC` 替换当前手写的远场和近场汇聚循环：

```cpp
ProcessMLFMM::lexpPEC(w,
    octreeNodes_[maxlvl],
    Vfmj,
    Bm,
    tN_k1,
    pN_k1,
    WGL_k1[0],
    WGL_phi_k1[0],
    maxlvl,
    k_,
    eta_,
    Z_near,
    Z_near_id,
    x);
```

这样可以保证 PEC 路径明确使用同一个 `lexpPEC` 权重入口。

### 6. PMCHWT 路径暂不混入本次改动

`MLFMM::applyAx` 的 PMCHWT 分支目前使用局部 `integrateL` 和 `integrateK` lambda，而不是 `lexpPEC`。本次目标限定为“权重只在 `lexpPEC` 中乘以 `lexpKpDot`”，因此建议先只改 PEC 路径。

如果后续也要修正 PMCHWT 的积分权重，应另行设计 `integrateL/integrateK` 的权重入口，避免和本次 PEC 改动混在一起。

## 验证方案

1. 全仓搜索确认 `m2l` 不再接收或使用 `WGL`、`WGL_phi`、`phi_level`。
2. 全仓搜索确认 `WGL[tt] * WGL_phi` 只出现在 `lexpPEC` 的积分累加处。
3. 编译 `MLFMA.sln`，确认所有 `m2l` 调用点和 `lexpPEC` 签名已同步。
4. 用现有输入运行一次 PEC 模式，对比修改前后输出量级，确认没有出现明显的重复乘权重或漏乘权重现象。

## 预期影响

- `Bm` 将只保存未经积分权重加权的局部展开系数。
- 权重位置更接近数学积分定义，集中在最终角谱求和处。
- `m2l` 接口更干净，减少无效参数和误用风险。
- PEC 路径会真正复用 `ProcessMLFMM::lexpPEC`，后续维护点更少。
