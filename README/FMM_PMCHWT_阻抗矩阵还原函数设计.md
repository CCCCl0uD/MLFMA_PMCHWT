# FMM/MLFMM PMCHWT 阻抗矩阵还原函数设计

## 目标

当前 MoM 的 PMCHWT 会直接生成完整的 `2N x 2N` 阻抗矩阵 `Z_mom`，而 FMM/MLFMM 的 PMCHWT 在迭代求解中只隐式实现矩阵-向量乘：

```text
w = Z_fmm * x = Z_far * x + Z_near * x
```

为了和 MoM 的 PMCHWT 阻抗矩阵逐项对比，需要增加一个调试/验证用函数，把 FMM/MLFMM 中隐式的近场和远场贡献还原成完整矩阵。

推荐思路是：不要重新写一套双重 RWG 积分或远场公式，而是把现有 PMCHWT matrix-vector product 抽成可复用的 `apply` 函数，然后对每个单位基向量 `e_j` 计算一次 `Z * e_j`。得到的结果向量就是完整矩阵第 `j` 列。

## 推荐接口

建议新增一个独立工具头文件，例如：

```text
Code/ZReconstruct.h
```

接口设计如下：

```cpp
#pragma once

#include <complex>
#include <string>
#include <vector>

namespace ZReconstruct {

enum class Part {
    NearOnly,
    FarOnly,
    Full
};

struct Options {
    Part part = Part::Full;
    bool showProgress = true;
    bool exportTxt = false;
    std::string exportPath;
};

using DenseZ = std::vector<std::vector<std::complex<double>>>;

template <class Solver>
DenseZ reconstructPMCHWTDense(Solver& solver, const Options& options = {});

template <class Solver>
void exportPMCHWTDense(
    Solver& solver,
    const std::string& path,
    const Options& options = {});

} // namespace ZReconstruct
```

其中 `Solver` 可以是 `FMM` 或 `MLFMM`。函数要求 solver 暴露一个统一的 PMCHWT 矩阵-向量乘接口：

```cpp
void applyPMCHWT(
    const std::complex<double>* x,
    std::complex<double>* w,
    ZReconstruct::Part part);
```

## 函数行为

设 `N = solver.row`，PMCHWT 未知量总数为 `N2 = 2 * N`。

未知量排列必须和现有 FMM/MLFMM 求解器一致：

```text
x[0:N)     : 电流相关未知量
x[N:2N)   : 磁流相关未知量
w[0:N)     : 第一组测试方程
w[N:2N)   : 第二组测试方程
```

完整矩阵按列还原：

```cpp
template <class Solver>
ZReconstruct::DenseZ ZReconstruct::reconstructPMCHWTDense(
    Solver& solver,
    const Options& options)
{
    if (solver.integralEquType_ != 2) {
        throw std::runtime_error("reconstructPMCHWTDense requires PMCHWT.");
    }

    const int N = solver.row;
    const int N2 = 2 * N;

    DenseZ Z(N2, std::vector<std::complex<double>>(N2, {0.0, 0.0}));
    std::vector<std::complex<double>> x(N2, {0.0, 0.0});
    std::vector<std::complex<double>> w(N2, {0.0, 0.0});

    for (int col = 0; col < N2; ++col) {
        std::fill(x.begin(), x.end(), std::complex<double>(0.0, 0.0));
        std::fill(w.begin(), w.end(), std::complex<double>(0.0, 0.0));

        x[col] = std::complex<double>(1.0, 0.0);
        solver.applyPMCHWT(x.data(), w.data(), options.part);

        for (int row = 0; row < N2; ++row) {
            Z[row][col] = w[row];
        }
    }

    return Z;
}
```

这样还原出来的 `Z` 与 solver 实际迭代时使用的矩阵严格一致，包括远场近似、近场精确项、当前未知量缩放和符号约定。

## `applyPMCHWT` 拆分设计

现有代码中，FMM/MLFMM 的 PMCHWT matvec 写在 `matrix_solver()` 的局部 lambda `ax` 内部。建议把这部分拆成三个成员函数：

```cpp
void applyPMCHWTFar(
    const std::complex<double>* x,
    std::complex<double>* w) const;

void applyPMCHWTNear(
    const std::complex<double>* x,
    std::complex<double>* w) const;

void applyPMCHWT(
    const std::complex<double>* x,
    std::complex<double>* w,
    ZReconstruct::Part part) const;
```

`applyPMCHWT()` 负责清零输出并组合近远场：

```cpp
void FMM::applyPMCHWT(
    const std::complex<double>* x,
    std::complex<double>* w,
    ZReconstruct::Part part) const
{
    const int N2 = 2 * row;
    std::fill(w, w + N2, std::complex<double>(0.0, 0.0));

    if (part == ZReconstruct::Part::FarOnly ||
        part == ZReconstruct::Part::Full) {
        applyPMCHWTFar(x, w);
    }

    if (part == ZReconstruct::Part::NearOnly ||
        part == ZReconstruct::Part::Full) {
        applyPMCHWTNear(x, w);
    }
}
```

`applyPMCHWTNear()` 可以直接复用当前近场 CSR 数据：

```cpp
void FMM::applyPMCHWTNear(
    const std::complex<double>* x,
    std::complex<double>* w) const
{
    const int N2 = 2 * row;

#pragma omp parallel for schedule(dynamic)
    for (int rowIdx = 0; rowIdx < N2; ++rowIdx) {
        std::complex<double> sum(0.0, 0.0);
        for (int j = 0; j < static_cast<int>(Z_near[rowIdx].size()); ++j) {
            sum += x[Z_near_id[rowIdx][j]] * Z_near[rowIdx][j];
        }
        w[rowIdx] += sum;
    }
}
```

`applyPMCHWTFar()` 则从当前 `matrix_solver()` 的 PMCHWT 分支中搬出远场部分：

- FMM：搬出 `aggregate -> translate -> integrateL/integrateK` 两个介质空间 `k1/k2` 的部分。
- MLFMM：搬出四个 pass：`k1/J`、`k1/M`、`k2/J`、`k2/M`，即当前对 `mexp -> m2m -> m2l -> l2l -> lexpPMCHWT` 的调用。
- 注意：不要在 `applyPMCHWTFar()` 内加入 `Z_near`，近场只由 `applyPMCHWTNear()` 负责。

改完以后，原来的 `matrix_solver()` 中 `ax` 可以变成：

```cpp
auto ax = [this](int n, const std::complex<double>* x, std::complex<double>* w) {
    if (integralEquType_ == 2) {
        if (n != 2 * row) {
            throw std::invalid_argument("PMCHWT matrix size must be 2 * row.");
        }
        applyPMCHWT(x, w, ZReconstruct::Part::Full);
        return;
    }

    // 保留原 PEC EFIE/CFIE 分支
};
```

## 输出格式

建议导出为三列文本，便于和 MoM 的 `_Zmatrix.txt` 做脚本对比：

```text
row    col    real    imag
0      0      ...     ...
0      1      ...     ...
```

对应导出函数：

```cpp
template <class Solver>
void ZReconstruct::exportPMCHWTDense(
    Solver& solver,
    const std::string& path,
    const Options& options)
{
    DenseZ Z = reconstructPMCHWTDense(solver, options);

    std::ofstream out(path);
    out << "row\tcol\treal\timag\n";

    for (int i = 0; i < static_cast<int>(Z.size()); ++i) {
        for (int j = 0; j < static_cast<int>(Z[i].size()); ++j) {
            out << i << '\t'
                << j << '\t'
                << Z[i][j].real() << '\t'
                << Z[i][j].imag() << '\n';
        }
    }
}
```

如果矩阵较大，也可以增加稀疏/阈值导出：

```cpp
double dropTol = 0.0; // abs(Zij) <= dropTol 时不输出
```

但用于和 MoM 做逐项误差分析时，建议先导出完整 dense 矩阵。

## 对比建议

建议至少导出三份矩阵：

```text
FMM_PMCHWT_Z_near.txt
FMM_PMCHWT_Z_far.txt
FMM_PMCHWT_Z_full.txt
```

其中：

- `NearOnly` 应该和 MoM 中同一近邻集合对应的直接积分项一致。
- `FarOnly` 是 FMM/MLFMM 近似远场项，通常不会和 MoM 逐项完全一致，但误差应随多极展开阶数、角向积分点、树划分参数改善。
- `Full` 是最终用于迭代求解的等效阻抗矩阵，最适合与 MoM 的完整 `Z_mom` 做整体误差统计。

推荐误差指标：

```text
maxAbsError = max_ij |Z_fmm(i,j) - Z_mom(i,j)|
relFroError = ||Z_fmm - Z_mom||_F / ||Z_mom||_F
blockRelError11 = ||Z11_fmm - Z11_mom||_F / ||Z11_mom||_F
blockRelError12 = ||Z12_fmm - Z12_mom||_F / ||Z12_mom||_F
blockRelError21 = ||Z21_fmm - Z21_mom||_F / ||Z21_mom||_F
blockRelError22 = ||Z22_fmm - Z22_mom||_F / ||Z22_mom||_F
```

## 注意事项

1. 这个函数主要用于小规模算例调试。完整还原 `2N x 2N` 矩阵需要做 `2N` 次 matvec，内存也是 `O((2N)^2)`。
2. 对比 MoM 前必须确认未知量缩放一致。当前 FMM 代码注释中提到 PMCHWT 的电流未知量可能使用 `eta1` 缩放；若 MoM 的 `Z_mom` 未采用同样缩放，需要先统一变量定义，否则矩阵块会出现系统性比例差异。
3. `MoM::exportZ()` 当前按 `row` 输出行数；PMCHWT 的完整矩阵是 `2 * row` 行。如果要直接导出完整 MoM PMCHWT 矩阵，建议把循环行数改为 `Z_mom.size()`。
4. `FMM` 和 `MLFMM` 的 `applyPMCHWTFar()` 应尽量只搬移现有 `matrix_solver()` 中已经验证过的 PMCHWT 远场逻辑，不要重写公式。这样还原矩阵才代表实际迭代器正在使用的矩阵。

