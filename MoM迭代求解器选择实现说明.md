# 在 MoM 中通过 `selectMatrixSolver` 选择 GMRES 或 CGS

## 1. 目标与当前代码结构

`main.cpp` 中已经定义：

```cpp
int selectMatrixSolver = 0; // 0 ==> GMRES; 1 ==> CGS
```

但这个变量目前没有传给 `MoM`，所以它不能影响实际求解过程。

当前 MoM 的调用链是：

```text
main.cpp
  └─ 构造 MoM
       └─ mom_Dual_* / mom_Mono_*
            └─ RCSUtils::compute*（RCS.h）
                 └─ solver.mgmres_solver(...)
                      └─ MGMRES::mgmres(...)
```

最小且稳妥的实现方式是：

1. 把 `selectMatrixSolver` 传入 `MoM` 构造函数；
2. 在 `MoM` 中保存该值；
3. 保留 `mgmres_solver(...)` 这个现有函数名，在函数内部选择 GMRES 或 CGS。

这里暂时保留 `mgmres_solver` 这个略显局限的旧名称，是为了不修改公共模板文件 `RCS.h`。`FMM` 和 `MLFMM` 也提供了同名函数，直接改名会扩大修改范围。

---

## 2. 修改 `main.cpp`

### 2.1 在构造 MoM 前校验选择值

建议在 `if (selectAlgorithm == 0)` 分支内、构造 `MoM` 之前加入：

```cpp
if (selectMatrixSolver != 0 && selectMatrixSolver != 1) {
    std::cerr << "Error: Invalid selectMatrixSolver value: "
              << selectMatrixSolver
              << " (0 = GMRES, 1 = CGS)" << std::endl;
    return 3;
}
```

这样错误值不会悄悄回退到某一种算法，也不会等到阻抗矩阵已经构造后才暴露问题。

### 2.2 把选择值传给 `MoM`

原代码：

```cpp
MoM mom(cfg, selectIntegralEqu, selectMono_Dual, polarization,
    octree.octreeNodes_, rwgBases, octree.maxLevel_, gausspoint,
    E0, wave);
```

改为：

```cpp
MoM mom(cfg, selectIntegralEqu, selectMatrixSolver,
    selectMono_Dual, polarization,
    octree.octreeNodes_, rwgBases, octree.maxLevel_, gausspoint,
    E0, wave);
```

修改后的完整 MoM 分支可以写成：

```cpp
if (selectAlgorithm == 0) {
    if (selectMatrixSolver != 0 && selectMatrixSolver != 1) {
        std::cerr << "Error: Invalid selectMatrixSolver value: "
                  << selectMatrixSolver
                  << " (0 = GMRES, 1 = CGS)" << std::endl;
        return 3;
    }

    auto startMoM = std::chrono::high_resolution_clock::now();
    MoM mom(cfg, selectIntegralEqu, selectMatrixSolver,
        selectMono_Dual, polarization,
        octree.octreeNodes_, rwgBases, octree.maxLevel_, gausspoint,
        E0, wave);
    auto endMoM = std::chrono::high_resolution_clock::now();
    auto durationMoM = std::chrono::duration_cast<std::chrono::microseconds>(
        endMoM - startMoM);
    std::cout << "<<<< Info: MoM done, elapsed: "
              << durationMoM.count() / 1000.0
              << " ms >>>> (main.cpp)\n" << std::endl;
}
```

---

## 3. 修改 `MoM.h`

### 3.1 修改构造函数声明

原声明：

```cpp
MoM(const RCSExportConfig& cfg, const int selectIntegralEqu,
    const std::string selectMono_Dual, const std::string pol_wave,
    const std::vector<std::vector<OCTree::Node*>>& octreeNodes,
    const std::vector<RWGBase>& rwgs, const int maxLevel_,
    const gaussPoints& gausspoint, const double E0,
    const EMSource& wave);
```

改为：

```cpp
MoM(const RCSExportConfig& cfg, const int selectIntegralEqu,
    const int selectMatrixSolver,
    const std::string selectMono_Dual, const std::string pol_wave,
    const std::vector<std::vector<OCTree::Node*>>& octreeNodes,
    const std::vector<RWGBase>& rwgs, const int maxLevel_,
    const gaussPoints& gausspoint, const double E0,
    const EMSource& wave);
```

### 3.2 增加成员变量

在现有的 `integralEquType_` 后加入：

```cpp
const int integralEquType_;
const int matrixSolverType_; // 0 = GMRES, 1 = CGS
```

### 3.3 求解函数声明保持不变

仍然保留：

```cpp
void mgmres_solver(int n, std::complex<double> x[],
    std::complex<double> rhs[], int itr_max, int mr,
    double tol_abs, double tol_rel);
```

虽然它会同时负责 GMRES 和 CGS，但保留名称可以让 `RCS.h` 的四处求解调用以及 FMM/MLFMA 接口完全不变。

---

## 4. 修改 `MoM.cpp`

### 4.1 引入 CGS 头文件

在文件头部加入：

```cpp
#include "CGS.h"
```

头部最终应类似：

```cpp
#include "Zmartix.h"
#include "Vm.h"
#include "GMRES.h"
#include "CGS.h"
#include "RCS.h"
#include "MoM.h"
```

### 4.2 修改构造函数定义及初始化列表

原定义的前半部分：

```cpp
MoM::MoM(
    const RCSExportConfig& cfg, const int selectIntegralEqu,
    const std::string selectMono_Dual, const std::string pol_wave,
    const std::vector<std::vector<OCTree::Node*>>& octreeNodes,
    const std::vector<RWGBase>& rwgs, const int maxLevel_,
    const gaussPoints& gausspoint, const double E0,
    const EMSource& wave)
    : octreeNodes_(octreeNodes), rwgs(rwgs), maxLevel_(maxLevel_),
      gausspoint(gausspoint), E0(E0), wave(wave),
      integralEquType_(selectIntegralEqu)
```

改为：

```cpp
MoM::MoM(
    const RCSExportConfig& cfg, const int selectIntegralEqu,
    const int selectMatrixSolver,
    const std::string selectMono_Dual, const std::string pol_wave,
    const std::vector<std::vector<OCTree::Node*>>& octreeNodes,
    const std::vector<RWGBase>& rwgs, const int maxLevel_,
    const gaussPoints& gausspoint, const double E0,
    const EMSource& wave)
    : octreeNodes_(octreeNodes),
      rwgs(rwgs),
      maxLevel_(maxLevel_),
      E0(E0),
      gausspoint(gausspoint),
      wave(wave),
      integralEquType_(selectIntegralEqu),
      matrixSolverType_(selectMatrixSolver)
```

初始化列表的书写顺序与 `MoM.h` 中成员的声明顺序一致，可避免 MSVC 的成员初始化顺序警告。

### 4.3 用分派逻辑替换 `mgmres_solver` 的函数体

用下面的完整实现替换当前 `MoM::mgmres_solver(...)`：

```cpp
void MoM::mgmres_solver(
    int n, std::complex<double> x[], std::complex<double> rhs[],
    int itr_max, int mr, double tol_abs, double tol_rel)
{
    // GMRES 和 CGS 共用同一个 MoM 矩阵-向量乘算子 w = Z_mom * x。
    auto ax = [this](
        int n,
        const std::complex<double>* x,
        std::complex<double>* w) {
#pragma omp parallel for
        for (int i = 0; i < n; ++i) {
            std::complex<double> sum = 0.0;
            for (int j = 0; j < n; ++j) {
                sum += Z_mom[i][j] * x[j];
            }
            w[i] = sum;
        }
    };

    bool converged = false;

    switch (matrixSolverType_) {
    case 0:
        std::cout << "Matrix solver: GMRES\n";
        converged = MGMRES::mgmres(
            n, x, rhs, itr_max, mr, tol_abs, tol_rel, ax);
        if (!converged) {
            std::cerr << "GMRES did not converge.\n";
        }
        break;

    case 1:
        std::cout << "Matrix solver: CGS\n";
        converged = MCGS::CGS::solve(
            n, x, rhs, itr_max, tol_rel, ax);
        if (!converged) {
            std::cerr << "CGS did not converge or encountered breakdown.\n";
        }
        break;

    default:
        throw std::invalid_argument(
            "Invalid matrixSolverType_ (must be 0 for GMRES or 1 for CGS)");
    }
}
```

`ax` 没有改变，两个求解器实际求解的都是：

```text
Z_mom · x = rhs
```

因此 EFIE、CFIE 和 PMCHWT 不需要分别编写选择逻辑；它们都会通过同一个入口自动选择求解器。PMCHWT 传入的 `n` 是 `2 * row`，而 `Z_mom` 也是对应的 `2N × 2N` 矩阵，该乘法逻辑仍然成立。

---

## 5. 两种求解器参数如何对应

现有 GMRES 接口是：

```cpp
MGMRES::mgmres(
    n, x, rhs,
    itr_max, mr,
    tol_abs, tol_rel,
    ax);
```

现有 CGS 接口是：

```cpp
MCGS::CGS::solve(
    n, x, rhs,
    max_iter, tol,
    ax);
```

参数对应关系如下：

| MoM 参数 | GMRES | CGS |
|---|---|---|
| `itr_max` | 最大总迭代次数 | 最大迭代次数 |
| `mr` | 重启时 Krylov 子空间大小 | 不使用 |
| `tol_abs` | 绝对残差阈值 | 当前实现不使用 |
| `tol_rel` | 相对残差阈值 | 作为 `tol` 使用 |

这是因为当前 `CGS.h` 的收敛判据只有：

```cpp
RSS = ||r|| / ||b||;
RSS < tol
```

所以调用 CGS 时使用 `tol_rel` 是语义一致的选择，不能简单使用 `tol_abs`。`mr` 是 GMRES 重启机制专用参数，对 CGS 没有意义。

---

## 6. 使用方式

选择 GMRES：

```cpp
int selectMatrixSolver = 0;
```

选择 CGS：

```cpp
int selectMatrixSolver = 1;
```

程序运行时会分别输出：

```text
Matrix solver: GMRES
```

或：

```text
Matrix solver: CGS
```

对于单站计算，每个入射角都会重新求解一次，因此这条提示也会出现多次，这是正常现象。

---

## 7. Visual Studio 项目中的可选整理

`CGS.h` 是头文件，只要 `MoM.cpp` 中能通过 `#include "CGS.h"` 找到它，就不必加入编译源文件列表，编译仍可成功。

不过当前 `MLFMA.vcxproj` 只列出了 `GMRES.h`，没有列出 `CGS.h`。如果希望它显示在 Visual Studio 的“头文件”节点中，可以进行以下可选修改。

在 `MLFMA.vcxproj` 的 `<ClInclude>` 列表中加入：

```xml
<ClInclude Include="CGS.h" />
```

在 `MLFMA.vcxproj.filters` 的头文件分组中加入：

```xml
<ClInclude Include="CGS.h">
  <Filter>头文件</Filter>
</ClInclude>
```

这两处只影响项目展示和文件管理，不影响上述选择功能本身。

---

## 8. 验证建议

完成代码修改后，建议依次进行以下验证：

1. 使用 `selectMatrixSolver = 0` 编译运行，确认输出 `Matrix solver: GMRES`，并与修改前的 MoM 结果对比；
2. 使用 `selectMatrixSolver = 1` 编译运行，确认输出 `Matrix solver: CGS`；
3. 分别测试 PEC/EFIE、PEC/CFIE 和介质/PMCHWT；
4. 检查最终相对残差是否低于 `tol_rel`；
5. 对同一模型比较两种算法得到的电流系数或 RCS，允许存在迭代误差量级内的小差异；
6. 将选择值改成 `2`，确认程序明确报错并退出。

CGS 对某些非厄米、病态的 MoM 矩阵可能发生 `rho` 为零或分母为零的算法击穿；这不是选择逻辑错误。当前 `CGS.h` 会返回 `false` 并输出 breakdown 信息。GMRES 通常更稳健，而 CGS 的内存需求较低，但收敛过程可能更不平滑。

---

## 9. 最终改动范围

实现该功能所需的必要源码改动只有：

```text
main.cpp  ：校验选择值，并将其传入 MoM
MoM.h     ：构造函数增加参数，类中增加选择成员
MoM.cpp   ：包含 CGS.h，并在统一求解入口中分派
```

不需要修改：

```text
RCS.h
GMRES.h
CGS.h
FMM.cpp / FMM.h
MLFMM.cpp / MLFMM.h
```

本文档仅给出实现代码和说明，未修改任何现有 C++ 源文件。
