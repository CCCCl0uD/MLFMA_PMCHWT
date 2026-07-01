# FMM/MLFMM PMCHWT 阻抗矩阵还原详细实现代码

## 实现目标

本文档给出一套可以直接落地到当前工程的实现方案，用于把 FMM/MLFMM 的 PMCHWT 隐式矩阵-向量乘还原成完整 `2N x 2N` 阻抗矩阵。

核心思想：

```text
Z(:, j) = Z * e_j
```

也就是对每一个单位基向量 `e_j` 调用一次现有 FMM/MLFMM 的 PMCHWT matvec，输出向量就是完整阻抗矩阵第 `j` 列。

建议新增或修改这些文件：

```text
Code/ZReconstruct.h
Code/FMM.h
Code/FMM.cpp
Code/MLFMM.h
Code/MLFMM.cpp
Code/MoM.cpp
```

## 1. 新增 `Code/ZReconstruct.h`

这个头文件提供矩阵还原和文本导出函数。它不关心具体 solver 是 `FMM` 还是 `MLFMM`，只要求 solver 暴露：

```cpp
void applyPMCHWT(
    const std::complex<double>* x,
    std::complex<double>* w,
    ZReconstruct::Part part);
```

完整代码如下：

```cpp
// ZReconstruct.h
#pragma once
#ifndef ZRECONSTRUCT_H
#define ZRECONSTRUCT_H

#include <algorithm>
#include <cmath>
#include <complex>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
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
    int progressStep = 1;
    double dropTol = 0.0;
};

using DenseZ = std::vector<std::vector<std::complex<double>>>;

inline const char* partName(Part part)
{
    switch (part) {
    case Part::NearOnly:
        return "near";
    case Part::FarOnly:
        return "far";
    case Part::Full:
        return "full";
    default:
        return "unknown";
    }
}

template <class Solver>
DenseZ reconstructPMCHWTDense(Solver& solver, const Options& options = {})
{
    if (solver.integralEquType_ != 2) {
        throw std::runtime_error("reconstructPMCHWTDense requires PMCHWT.");
    }

    const int N = solver.row;
    const int N2 = 2 * N;

    if (N <= 0) {
        throw std::runtime_error("reconstructPMCHWTDense got invalid solver.row.");
    }

    DenseZ Z(N2, std::vector<std::complex<double>>(N2, {0.0, 0.0}));
    std::vector<std::complex<double>> x(N2, {0.0, 0.0});
    std::vector<std::complex<double>> w(N2, {0.0, 0.0});

    const int step = std::max(1, options.progressStep);

    for (int col = 0; col < N2; ++col) {
        std::fill(x.begin(), x.end(), std::complex<double>(0.0, 0.0));
        std::fill(w.begin(), w.end(), std::complex<double>(0.0, 0.0));

        x[col] = std::complex<double>(1.0, 0.0);
        solver.applyPMCHWT(x.data(), w.data(), options.part);

        for (int row = 0; row < N2; ++row) {
            Z[row][col] = w[row];
        }

        if (options.showProgress && ((col + 1) % step == 0 || col + 1 == N2)) {
            std::cout << "Info: reconstruct PMCHWT Z_"
                << partName(options.part) << " column "
                << (col + 1) << " / " << N2 << std::endl;
        }
    }

    return Z;
}

template <class Solver>
void exportPMCHWTDense(
    Solver& solver,
    const std::string& path,
    const Options& options = {})
{
    DenseZ Z = reconstructPMCHWTDense(solver, options);

    std::ofstream out(path);
    if (!out.is_open()) {
        throw std::runtime_error("Failed to open PMCHWT Z export file: " + path);
    }

    out << std::setprecision(17);
    out << "row\tcol\treal\timag\n";

    for (int i = 0; i < static_cast<int>(Z.size()); ++i) {
        for (int j = 0; j < static_cast<int>(Z[i].size()); ++j) {
            const std::complex<double>& zij = Z[i][j];
            if (options.dropTol > 0.0 && std::abs(zij) <= options.dropTol) {
                continue;
            }

            out << i << '\t'
                << j << '\t'
                << zij.real() << '\t'
                << zij.imag() << '\n';
        }
    }

    std::cout << "Info: exported PMCHWT Z_"
        << partName(options.part) << " to " << path << std::endl;
}

} // namespace ZReconstruct

#endif // ZRECONSTRUCT_H
```

## 2. 修改 `Code/FMM.h`

在 `FMM.h` 中加入 `ZReconstruct.h`，并声明三个 PMCHWT matvec 函数。

建议把 include 加在已有 include 后：

```cpp
#include "ZReconstruct.h"
```

在 `public:` 区域加入：

```cpp
void applyPMCHWT(
    const std::complex<double>* x,
    std::complex<double>* w,
    ZReconstruct::Part part);

void applyPMCHWTNear(
    const std::complex<double>* x,
    std::complex<double>* w) const;

void applyPMCHWTFar(
    const std::complex<double>* x,
    std::complex<double>* w) const;
```

加入后的关键片段类似：

```cpp
#include "OCTree.h"
#include "RCSExportConfig.h"
#include "EMSource.h"
#include "GaussPoints.h"
#include "ZReconstruct.h"

class FMM {
public:
    FMM(...);

    // 省略已有成员

    size_t computeMem();
    void matrix_solver(int n, std::complex<double> x[], std::complex<double> rhs[],
        int itr_max, int mr, double tol_abs, double tol_rel);

    void applyPMCHWT(
        const std::complex<double>* x,
        std::complex<double>* w,
        ZReconstruct::Part part);

    void applyPMCHWTNear(
        const std::complex<double>* x,
        std::complex<double>* w) const;

    void applyPMCHWTFar(
        const std::complex<double>* x,
        std::complex<double>* w) const;

private:
    // 保留已有 private 成员
};
```

## 3. 在 `Code/FMM.cpp` 中实现 FMM 的还原用 matvec

把以下代码加到 `FMM.cpp` 中，建议放在 `computeMem()` 后、`matrix_solver()` 前。

```cpp
void FMM::applyPMCHWT(
    const std::complex<double>* x,
    std::complex<double>* w,
    ZReconstruct::Part part)
{
    if (integralEquType_ != 2) {
        throw std::runtime_error("FMM::applyPMCHWT requires PMCHWT.");
    }

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

void FMM::applyPMCHWTNear(
    const std::complex<double>* x,
    std::complex<double>* w) const
{
    const int N = row;
    const int N2 = 2 * N;

    if (static_cast<int>(Z_near.size()) != N2 ||
        static_cast<int>(Z_near_id.size()) != N2) {
        throw std::runtime_error("FMM PMCHWT near matrix must have 2 * row rows.");
    }

#pragma omp parallel for schedule(dynamic)
    for (int rowIdx = 0; rowIdx < N2; ++rowIdx) {
        std::complex<double> sum(0.0, 0.0);
        for (int j = 0; j < static_cast<int>(Z_near[rowIdx].size()); ++j) {
            sum += x[Z_near_id[rowIdx][j]] * Z_near[rowIdx][j];
        }
        w[rowIdx] += sum;
    }
}

void FMM::applyPMCHWTFar(
    const std::complex<double>* x,
    std::complex<double>* w) const
{
    const int N = row;
    const int nodenum = static_cast<int>(octreeNodes_[maxLevel_].size());

    const int K1 = static_cast<int>(kp_lvl_k1[0].size());
    const int K2 = static_cast<int>(kp_lvl_k2[0].size());

    const int tN1 = static_cast<int>(theta_level_k1[0].size());
    const int pN1 = static_cast<int>(phi_level_k1[0].size());
    const int tN2 = static_cast<int>(theta_level_k2[0].size());
    const int pN2 = static_cast<int>(phi_level_k2[0].size());

    const std::complex<double> farL1 = wave.k1() * wave.k1() / (16.0 * Pi * Pi);
    const std::complex<double> farL2 = wave.k2() * wave.k2() / (16.0 * Pi * Pi);
    const std::complex<double> etai = wave.etai();

    auto aggregate = [&](const std::vector<std::vector<Complex3D>>& Vsm,
        int xOffset, int K, std::vector<std::vector<Complex3D>>& Sm) {
#pragma omp parallel for schedule(dynamic)
            for (int nodeIdx = 0; nodeIdx < nodenum; ++nodeIdx) {
                OCTree::Node* node = octreeNodes_[maxLevel_][nodeIdx];
                const int rwgNum = static_cast<int>(node->rwgIndices.size());

                for (int r = 0; r < rwgNum; ++r) {
                    const int rwgid = node->rwgIndices[r]->rwgid;
                    const std::complex<double> coeff = x[xOffset + rwgid];

                    for (int kpi = 0; kpi < K; ++kpi) {
                        Sm[nodeIdx][kpi].x += Vsm[rwgid][kpi].x * coeff;
                        Sm[nodeIdx][kpi].y += Vsm[rwgid][kpi].y * coeff;
                        Sm[nodeIdx][kpi].z += Vsm[rwgid][kpi].z * coeff;
                    }
                }
            }
        };

    auto translate = [&](const std::vector<std::vector<std::complex<double>>>& TF_all,
        int K, const std::vector<std::vector<Complex3D>>& Sm,
        std::vector<std::vector<Complex3D>>& Bm) {
#pragma omp parallel for schedule(dynamic)
            for (int nodeIdx = 0; nodeIdx < nodenum; ++nodeIdx) {
                const int farnum = static_cast<int>(node_far[nodeIdx].size());

                for (int f = 0; f < farnum; ++f) {
                    const int farId = node_far_id[nodeIdx][f];
                    const auto& TF = TF_all[node_farVec_id[nodeIdx][f]];

                    for (int kpi = 0; kpi < K; ++kpi) {
                        Bm[nodeIdx][kpi].x += TF[kpi] * Sm[farId][kpi].x;
                        Bm[nodeIdx][kpi].y += TF[kpi] * Sm[farId][kpi].y;
                        Bm[nodeIdx][kpi].z += TF[kpi] * Sm[farId][kpi].z;
                    }
                }
            }
        };

    auto integrateL = [](const std::vector<Complex3D>& V,
        const std::vector<Complex3D>& B, const std::vector<double>& WGL,
        double WGL_phi, int tN, int pN) -> std::complex<double> {
            std::complex<double> val(0.0, 0.0);
            for (int ti = 0; ti < tN; ++ti) {
                for (int pj = 0; pj < pN; ++pj) {
                    const int idx = ti * pN + pj;
                    const double d2k = WGL[ti] * WGL_phi;
                    val += d2k * dot(V[idx], B[idx]);
                }
            }
            return val;
        };

    auto integrateK = [](const std::vector<Complex3D>& V,
        const std::vector<Complex3D>& B, const std::vector<kp_Point>& kp,
        const std::vector<double>& WGL, double WGL_phi,
        int tN, int pN) -> std::complex<double> {
            std::complex<double> val(0.0, 0.0);
            for (int ti = 0; ti < tN; ++ti) {
                for (int pj = 0; pj < pN; ++pj) {
                    const int idx = ti * pN + pj;
                    const double d2k = WGL[ti] * WGL_phi;

                    const double khat[3] = {
                        kp[idx].x,
                        kp[idx].y,
                        kp[idx].z
                    };

                    Complex3D kCrossV;
                    cross(kCrossV, khat, V[idx]);
                    val += d2k * dot(kCrossV, B[idx]);
                }
            }
            return val;
        };

    // k1 external space
    {
        std::vector<std::vector<Complex3D>> Sm_J1(
            nodenum, std::vector<Complex3D>(K1, Complex3D()));
        std::vector<std::vector<Complex3D>> Bm_J1(
            nodenum, std::vector<Complex3D>(K1, Complex3D()));
        std::vector<std::vector<Complex3D>> Sm_M1(
            nodenum, std::vector<Complex3D>(K1, Complex3D()));
        std::vector<std::vector<Complex3D>> Bm_M1(
            nodenum, std::vector<Complex3D>(K1, Complex3D()));

        aggregate(Vsmi, 0, K1, Sm_J1);
        aggregate(Vsmi, N, K1, Sm_M1);

        translate(TF_fmm, K1, Sm_J1, Bm_J1);
        translate(TF_fmm, K1, Sm_M1, Bm_M1);

#pragma omp parallel for schedule(dynamic)
        for (int nodeIdx = 0; nodeIdx < nodenum; ++nodeIdx) {
            OCTree::Node* node = octreeNodes_[maxLevel_][nodeIdx];
            const int rwgNum = static_cast<int>(node->rwgIndices.size());

            for (int r = 0; r < rwgNum; ++r) {
                const int rwgid = node->rwgIndices[r]->rwgid;
                const auto& V1 = Vfmj[rwgid];

                const auto L1J = integrateL(
                    V1, Bm_J1[nodeIdx], WGL_k1[0], WGL_phi_k1[0], tN1, pN1);
                const auto K1J = integrateK(
                    V1, Bm_J1[nodeIdx], kp_lvl_k1[0], WGL_k1[0], WGL_phi_k1[0], tN1, pN1);
                const auto K1M = integrateK(
                    V1, Bm_M1[nodeIdx], kp_lvl_k1[0], WGL_k1[0], WGL_phi_k1[0], tN1, pN1);
                const auto L1M = integrateL(
                    V1, Bm_M1[nodeIdx], WGL_k1[0], WGL_phi_k1[0], tN1, pN1);

                w[rwgid] += farL1 * L1J;
                w[N + rwgid] += farL1 * K1J;
                w[rwgid] += -farL1 * K1M;
                w[N + rwgid] += farL1 * L1M;
            }
        }
    }

    // k2 internal space
    {
        std::vector<std::vector<Complex3D>> Sm_J2(
            nodenum, std::vector<Complex3D>(K2, Complex3D()));
        std::vector<std::vector<Complex3D>> Bm_J2(
            nodenum, std::vector<Complex3D>(K2, Complex3D()));
        std::vector<std::vector<Complex3D>> Sm_M2(
            nodenum, std::vector<Complex3D>(K2, Complex3D()));
        std::vector<std::vector<Complex3D>> Bm_M2(
            nodenum, std::vector<Complex3D>(K2, Complex3D()));

        aggregate(Vsmi2, 0, K2, Sm_J2);
        aggregate(Vsmi2, N, K2, Sm_M2);

        translate(TF_fmm2, K2, Sm_J2, Bm_J2);
        translate(TF_fmm2, K2, Sm_M2, Bm_M2);

#pragma omp parallel for schedule(dynamic)
        for (int nodeIdx = 0; nodeIdx < nodenum; ++nodeIdx) {
            OCTree::Node* node = octreeNodes_[maxLevel_][nodeIdx];
            const int rwgNum = static_cast<int>(node->rwgIndices.size());

            for (int r = 0; r < rwgNum; ++r) {
                const int rwgid = node->rwgIndices[r]->rwgid;
                const auto& V2 = Vfmj2[rwgid];

                const auto L2J = integrateL(
                    V2, Bm_J2[nodeIdx], WGL_k2[0], WGL_phi_k2[0], tN2, pN2);
                const auto K2J = integrateK(
                    V2, Bm_J2[nodeIdx], kp_lvl_k2[0], WGL_k2[0], WGL_phi_k2[0], tN2, pN2);
                const auto K2M = integrateK(
                    V2, Bm_M2[nodeIdx], kp_lvl_k2[0], WGL_k2[0], WGL_phi_k2[0], tN2, pN2);
                const auto L2M = integrateL(
                    V2, Bm_M2[nodeIdx], WGL_k2[0], WGL_phi_k2[0], tN2, pN2);

                w[rwgid] += etai * farL2 * L2J;
                w[N + rwgid] += farL2 * K2J;
                w[rwgid] += -farL2 * K2M;
                w[N + rwgid] += (farL2 / etai) * L2M;
            }
        }
    }
}
```

## 4. 替换 `FMM::matrix_solver()` 中 PMCHWT 分支

原来 `matrix_solver()` 的 PMCHWT 分支很长。完成上面的拆分后，可以把 `ax` 里 PMCHWT 的部分替换为：

```cpp
auto ax = [this](int n, const std::complex<double>* x, std::complex<double>* w) {
    for (int i = 0; i < n; ++i) {
        w[i] = std::complex<double>(0.0, 0.0);
    }

    if (integralEquType_ == 2) {
        const int N = row;
        if (n != 2 * N) {
            throw std::invalid_argument("FMM PMCHWT matrix_solver requires n == 2 * row.");
        }

        applyPMCHWT(x, w, ZReconstruct::Part::Full);
        return;
    }

    // 这里保留原来的 PEC EFIE/CFIE matvec 分支
};
```

注意：`applyPMCHWT()` 内部已经会清零 `w`，所以外层清零保留或删除都可以。为了和原代码风格一致，可以保留。

## 5. 修改 `Code/MLFMM.h`

和 `FMM.h` 类似，先增加：

```cpp
#include "ZReconstruct.h"
```

然后在 `public:` 中加入：

```cpp
void applyPMCHWT(
    const std::complex<double>* x,
    std::complex<double>* w,
    ZReconstruct::Part part);

void applyPMCHWTNear(
    const std::complex<double>* x,
    std::complex<double>* w) const;

void applyPMCHWTFar(
    const std::complex<double>* x,
    std::complex<double>* w);
```

`MLFMM::applyPMCHWTFar()` 不建议声明为 `const`，因为当前 MLFMM 远场 matvec 会复用并清空成员变量 `Sm/Bm/Sm2/Bm2`。

## 6. 在 `Code/MLFMM.cpp` 中实现 MLFMM 的还原用 matvec

建议放在 `computeMem()` 后、`matrix_solver()` 前。

```cpp
void MLFMM::applyPMCHWT(
    const std::complex<double>* x,
    std::complex<double>* w,
    ZReconstruct::Part part)
{
    if (integralEquType_ != 2) {
        throw std::runtime_error("MLFMM::applyPMCHWT requires PMCHWT.");
    }

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

void MLFMM::applyPMCHWTNear(
    const std::complex<double>* x,
    std::complex<double>* w) const
{
    const int N = row;
    const int N2 = 2 * N;

    if (static_cast<int>(Z_near.size()) != N2 ||
        static_cast<int>(Z_near_id.size()) != N2) {
        throw std::runtime_error("MLFMM PMCHWT near matrix must have 2 * row rows.");
    }

    const int finestN = static_cast<int>(octreeNodes_[maxLevel_].size());

#pragma omp parallel for schedule(dynamic)
    for (int nodeIdx = 0; nodeIdx < finestN; ++nodeIdx) {
        OCTree::Node* nd = octreeNodes_[maxLevel_][nodeIdx];
        for (int ri = 0; ri < static_cast<int>(nd->rwgIndices.size()); ++ri) {
            const int rwgid = nd->rwgIndices[ri]->rwgid;

            std::complex<double> sumJ(0.0, 0.0);
            for (int j = 0; j < static_cast<int>(Z_near[rwgid].size()); ++j) {
                sumJ += x[Z_near_id[rwgid][j]] * Z_near[rwgid][j];
            }
            w[rwgid] += sumJ;

            std::complex<double> sumM(0.0, 0.0);
            for (int j = 0; j < static_cast<int>(Z_near[N + rwgid].size()); ++j) {
                sumM += x[Z_near_id[N + rwgid][j]] * Z_near[N + rwgid][j];
            }
            w[N + rwgid] += sumM;
        }
    }
}

void MLFMM::applyPMCHWTFar(
    const std::complex<double>* x,
    std::complex<double>* w)
{
    const int N = row;
    const int maxlvl = maxLevel_;

    const int tN_k1 = static_cast<int>(theta_level_k1[0].size());
    const int pN_k1 = static_cast<int>(phi_level_k1[0].size());
    const int tN_k2 = static_cast<int>(theta_level_k2[0].size());
    const int pN_k2 = static_cast<int>(phi_level_k2[0].size());

    const std::complex<double> const1 =
        wave.k1() * wave.k1() * wave.eta1() / (16.0 * Pi * Pi);
    const std::complex<double> const2 =
        wave.k2() * wave.k2() / (16.0 * Pi * Pi);

    // Pass 1: k1 space, J source -> L1(J), K1(J)
    ProcessMLFMM::clearSmBm(Sm, Bm);
    ProcessMLFMM::mexp(Sm, octreeNodes_[maxLevel_], Vsmi,
        kp_lvl_k1[0], x, 0, maxlvl);
    ProcessMLFMM::m2m(Sm, octreeNodes_, kp_lvl_k1,
        interpol_k1, maxlvl, levelSpan, wave.k1());
    ProcessMLFMM::m2l(Sm, Bm, octreeNodes_, TF, kp_lvl_k1,
        WGL_k1, WGL_phi_k1, phi_level_k1, maxlvl, levelSpan);
    ProcessMLFMM::l2l(Bm, octreeNodes_, kp_lvl_k1,
        phi_level_k1, interpol_k1, maxlvl, levelSpan, wave.k1());
    ProcessMLFMM::lexpPMCHWT(w, const1, -const1, 0, N, N,
        octreeNodes_[maxLevel_], Vfmj, Bm, kp_lvl_k1[0],
        tN_k1, pN_k1, maxlvl);

    // Pass 2: k1 space, M source -> K1(M), L1(M)
    ProcessMLFMM::clearSmBm(Sm, Bm);
    ProcessMLFMM::mexp(Sm, octreeNodes_[maxLevel_], Vsmi,
        kp_lvl_k1[0], x, N, maxlvl);
    ProcessMLFMM::m2m(Sm, octreeNodes_, kp_lvl_k1,
        interpol_k1, maxlvl, levelSpan, wave.k1());
    ProcessMLFMM::m2l(Sm, Bm, octreeNodes_, TF, kp_lvl_k1,
        WGL_k1, WGL_phi_k1, phi_level_k1, maxlvl, levelSpan);
    ProcessMLFMM::l2l(Bm, octreeNodes_, kp_lvl_k1,
        phi_level_k1, interpol_k1, maxlvl, levelSpan, wave.k1());
    ProcessMLFMM::lexpPMCHWT(w, const1, const1, N, 0, N,
        octreeNodes_[maxLevel_], Vfmj, Bm, kp_lvl_k1[0],
        tN_k1, pN_k1, maxlvl);

    // Pass 3: k2 space, J source -> L2(J), K2(J)
    ProcessMLFMM::clearSmBm(Sm2, Bm2);
    ProcessMLFMM::mexp(Sm2, octreeNodes_[maxLevel_], Vsmi2,
        kp_lvl_k2[0], x, 0, maxlvl);
    ProcessMLFMM::m2m(Sm2, octreeNodes_, kp_lvl_k2,
        interpol_k2, maxlvl, levelSpan, wave.k2());
    ProcessMLFMM::m2l(Sm2, Bm2, octreeNodes_, TF2, kp_lvl_k2,
        WGL_k2, WGL_phi_k2, phi_level_k2, maxlvl, levelSpan);
    ProcessMLFMM::l2l(Bm2, octreeNodes_, kp_lvl_k2,
        phi_level_k2, interpol_k2, maxlvl, levelSpan, wave.k2());
    ProcessMLFMM::lexpPMCHWT(w,
        const2 * wave.eta2(), -const2 * wave.eta1(), 0, N, N,
        octreeNodes_[maxLevel_], Vfmj2, Bm2, kp_lvl_k2[0],
        tN_k2, pN_k2, maxlvl);

    // Pass 4: k2 space, M source -> K2(M), L2(M)
    ProcessMLFMM::clearSmBm(Sm2, Bm2);
    ProcessMLFMM::mexp(Sm2, octreeNodes_[maxLevel_], Vsmi2,
        kp_lvl_k2[0], x, N, maxlvl);
    ProcessMLFMM::m2m(Sm2, octreeNodes_, kp_lvl_k2,
        interpol_k2, maxlvl, levelSpan, wave.k2());
    ProcessMLFMM::m2l(Sm2, Bm2, octreeNodes_, TF2, kp_lvl_k2,
        WGL_k2, WGL_phi_k2, phi_level_k2, maxlvl, levelSpan);
    ProcessMLFMM::l2l(Bm2, octreeNodes_, kp_lvl_k2,
        phi_level_k2, interpol_k2, maxlvl, levelSpan, wave.k2());
    ProcessMLFMM::lexpPMCHWT(w,
        const2 * wave.eta2() * wave.epsilonR(), const2 * wave.eta1(),
        N, 0, N,
        octreeNodes_[maxLevel_], Vfmj2, Bm2, kp_lvl_k2[0],
        tN_k2, pN_k2, maxlvl);
}
```

重要检查点：上面 `l2l()` 的实参已经按当前工程中的 `ProcessMLFMM::l2l()` 声明写成带 `levelSpan` 的版本。你当前 `MLFMM.cpp` 原始调用是：

```cpp
ProcessMLFMM::l2l(Bm, octreeNodes_, kp_lvl_k1,
    phi_level_k1, interpol_k1, maxlvl, levelSpan, wave.k1());
```

也就是说，搬移代码时以当前工程中 `matrix_solver()` 里已经能编译的调用为准。

## 7. 替换 `MLFMM::matrix_solver()` 中 PMCHWT 分支

完成拆分后，`matrix_solver()` 中的 `ax` 可以改成：

```cpp
const int tN_k1 = static_cast<int>(theta_level_k1[0].size());
const int pN_k1 = static_cast<int>(phi_level_k1[0].size());
const int maxlvl = maxLevel_;

auto ax = [this, tN_k1, pN_k1, maxlvl](
    int n,
    const std::complex<double>* x,
    std::complex<double>* w) {

    for (int i = 0; i < n; ++i) {
        w[i] = std::complex<double>(0.0, 0.0);
    }

    if (integralEquType_ == 2) {
        if (n != 2 * row) {
            throw std::invalid_argument("MLFMM PMCHWT matrix_solver requires n == 2 * row.");
        }
        applyPMCHWT(x, w, ZReconstruct::Part::Full);
        return;
    }

    // 这里保留原来的 PEC EFIE/CFIE 分支
    ProcessMLFMM::clearSmBm(Sm, Bm);
    ProcessMLFMM::mexp(Sm, octreeNodes_[maxLevel_], Vsmi,
        kp_lvl_k1[0], x, 0, maxlvl);
    ProcessMLFMM::m2m(Sm, octreeNodes_, kp_lvl_k1,
        interpol_k1, maxlvl, levelSpan, wave.k1());
    ProcessMLFMM::m2l(Sm, Bm, octreeNodes_, TF, kp_lvl_k1,
        WGL_k1, WGL_phi_k1, phi_level_k1, maxlvl, levelSpan);
    ProcessMLFMM::l2l(Bm, octreeNodes_, kp_lvl_k1,
        phi_level_k1, interpol_k1, maxlvl, levelSpan, wave.k1());
    ProcessMLFMM::lexpPEC(w, octreeNodes_[maxLevel_], Vfmj, Bm,
        tN_k1, pN_k1, maxlvl,
        wave.k1(), wave.eta1(), Z_near, Z_near_id, x);
};
```

## 8. 在主流程中导出 FMM/MLFMM PMCHWT 矩阵

在需要调试的地方 include：

```cpp
#include "ZReconstruct.h"
```

如果要在 `main.cpp` 中构造 `FMM fmm(...)` 后立刻导出：

```cpp
if (selectIntegralEqu == 2) {
    ZReconstruct::Options opt;
    opt.showProgress = true;
    opt.progressStep = 10;

    opt.part = ZReconstruct::Part::NearOnly;
    ZReconstruct::exportPMCHWTDense(
        fmm,
        "FMM_PMCHWT_Z_near.txt",
        opt);

    opt.part = ZReconstruct::Part::FarOnly;
    ZReconstruct::exportPMCHWTDense(
        fmm,
        "FMM_PMCHWT_Z_far.txt",
        opt);

    opt.part = ZReconstruct::Part::Full;
    ZReconstruct::exportPMCHWTDense(
        fmm,
        "FMM_PMCHWT_Z_full.txt",
        opt);
}
```

如果是 `MLFMM mlfmm(...)`，调用方式完全一样：

```cpp
if (selectIntegralEqu == 2) {
    ZReconstruct::Options opt;
    opt.showProgress = true;
    opt.progressStep = 10;

    opt.part = ZReconstruct::Part::Full;
    ZReconstruct::exportPMCHWTDense(
        mlfmm,
        "MLFMM_PMCHWT_Z_full.txt",
        opt);
}
```

## 9. 修正 MoM 的 PMCHWT 完整矩阵导出

当前 `MoM::exportZ()` 使用：

```cpp
for (int i = 0; i < row; i++)
```

但 PMCHWT 的 `Z_mom` 是 `2 * row` 行。建议改为：

```cpp
void MoM::exportZ(const RCSExportConfig& cfg) {
    auto zFile = cfg.openFile("_Zmatrix", "txt");
    if (!zFile.is_open()) {
        std::cerr << "Error: Failed to open file for Z matrix export." << std::endl;
        return;
    }

    zFile << std::setprecision(17);
    zFile << "row\tcol\treal\timag\n";

    for (int i = 0; i < static_cast<int>(Z_mom.size()); ++i) {
        for (int j = 0; j < static_cast<int>(Z_mom[i].size()); ++j) {
            zFile << i << '\t'
                << j << '\t'
                << Z_mom[i][j].real() << '\t'
                << Z_mom[i][j].imag() << '\n';
        }
    }

    zFile.close();
    std::cout << "Z matrix export completed." << std::endl;
}
```

如果使用 `std::setprecision`，需要在 `MoM.cpp` 顶部加入：

```cpp
#include <iomanip>
```

## 10. 误差对比脚本思路

导出 `MoM_PMCHWT_Z_full.txt` 和 `FMM_PMCHWT_Z_full.txt` 后，可以用任意脚本计算误差。C++ 中也可以新增一个简单读取器，不过更建议先用 Python 或 MATLAB 做验证。

建议统计：

```text
maxAbsError = max_ij |Z_fmm(i,j) - Z_mom(i,j)|
relFroError = ||Z_fmm - Z_mom||_F / ||Z_mom||_F
Z11RelError = ||Z11_fmm - Z11_mom||_F / ||Z11_mom||_F
Z12RelError = ||Z12_fmm - Z12_mom||_F / ||Z12_mom||_F
Z21RelError = ||Z21_fmm - Z21_mom||_F / ||Z21_mom||_F
Z22RelError = ||Z22_fmm - Z22_mom||_F / ||Z22_mom||_F
```

PMCHWT 四个块的索引：

```text
Z11 = Z[0:N,     0:N]
Z12 = Z[0:N,     N:2N]
Z21 = Z[N:2N,    0:N]
Z22 = Z[N:2N,    N:2N]
```

## 11. 实现注意事项

1. 完整还原矩阵只适合小规模算例。若 RWG 数为 `N`，需要做 `2N` 次 matvec，并保存 `(2N)^2` 个复数。
2. `NearOnly` 输出的是当前 FMM/MLFMM 近邻判定下的精确积分块，不是 MoM 全矩阵。
3. `FarOnly` 输出的是 FMM/MLFMM 近似远场贡献，不应期望每一项都和 MoM 逐项完全一致。
4. `Full` 输出的是迭代求解器实际使用的等效矩阵，最适合与 MoM 完整矩阵做总体误差分析。
5. 如果 FMM 和 MoM 的 PMCHWT 未知量缩放不同，例如 FMM 中电流未知量使用 `eta1 * J`，则必须先统一变量定义，否则四个矩阵块会出现系统性比例差异。
6. 搬移 MLFMM 远场代码时，`mexp/m2m/m2l/l2l/lexpPMCHWT` 的参数顺序必须以当前工程能编译的原始调用为准。
