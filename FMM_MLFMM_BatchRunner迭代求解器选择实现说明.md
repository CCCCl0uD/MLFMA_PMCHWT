# FMM、MLFMM 与 BatchRunner 的 GMRES/CGS 选择实现

## 1. 实现目标

你已经在 `main.cpp` 和 `MoM` 中完成了以下设计：

```cpp
int selectMatrixSolver = 0; // 0 = GMRES; 1 = CGS
```

并把三种算法供 `RCS.h` 调用的统一接口命名为：

```cpp
void matrix_solver(
    int n,
    std::complex<double> x[],
    std::complex<double> rhs[],
    int itr_max,
    int mr,
    double tol_abs,
    double tol_rel);
```

接下来应让控制值按下面的路径传递：

```text
单次运行：main.cpp
  ├─ MoM
  ├─ FMM
  └─ MLFMM

批量运行：batch_config.txt
  └─ BatchConfig
       └─ SimTask
            └─ BatchRunner::executeSingleTask
                 ├─ MoM
                 ├─ FMM
                 └─ MLFMM
```

FMM 和 MLFMM 已经分别把完整矩阵作用封装在 `matrix_solver()` 内的 `ax` lambda 中。因此不要复制或修改庞大的近场/远场计算代码，只需保存选择值，并在 `ax` 定义之后分派 GMRES 或 CGS。

---

## 2. 修改 `FMM.h`

### 2.1 构造函数增加参数

将原声明：

```cpp
FMM(const RCSExportConfig& cfg, const int selectIntegralEqu,
    const std::string selectMono_Dual, const std::string pol_wave,
    const std::vector<std::vector<OCTree::Node*>>& octreeNodes,
    const std::vector<RWGBase>& rwgs, const int maxLevel_,
    const gaussPoints& gausspoint, const double E0,
    const EMSource& wave);
```

改为：

```cpp
FMM(const RCSExportConfig& cfg, const int selectIntegralEqu,
    const int selectMatrixSolver,
    const std::string selectMono_Dual, const std::string pol_wave,
    const std::vector<std::vector<OCTree::Node*>>& octreeNodes,
    const std::vector<RWGBase>& rwgs, const int maxLevel_,
    const gaussPoints& gausspoint, const double E0,
    const EMSource& wave);
```

### 2.2 保存选择值

在 `integralEquType_` 后加入：

```cpp
const int integralEquType_;
const int matrixSolverType_; // 0 = GMRES, 1 = CGS
```

现有 `matrix_solver(...)` 声明保持不变。

---

## 3. 修改 `FMM.cpp`

### 3.1 引入 CGS

在 `GMRES.h` 后加入：

```cpp
#include <stdexcept>

#include "GMRES.h"
#include "CGS.h"
```

### 3.2 修改构造函数定义

将构造函数头部和初始化列表改为：

```cpp
FMM::FMM(
    const RCSExportConfig& cfg,
    const int selectIntegralEqu,
    const int selectMatrixSolver,
    const std::string selectMono_Dual,
    const std::string pol_wave,
    const std::vector<std::vector<OCTree::Node*>>& octreeNodes,
    const std::vector<RWGBase>& rwgs,
    const int maxLevel_,
    const gaussPoints& gausspoint,
    const double E0,
    const EMSource& wave)
    : octreeNodes_(octreeNodes),
      rwgs(rwgs),
      maxLevel_(maxLevel_),
      E0(E0),
      gausspoint(gausspoint),
      wave(wave),
      integralEquType_(selectIntegralEqu),
      matrixSolverType_(selectMatrixSolver)
{
    // 原有构造函数主体保持不变
```

### 3.3 替换 `matrix_solver()` 末尾的 GMRES 固定调用

保留当前函数中从开头到下面这一行为止的全部代码：

```cpp
}; // end ax lambda
```

仅把其后的原代码：

```cpp
bool converged = MGMRES::mgmres(
    n, x, rhs, itr_max, mr, tol_abs, tol_rel, ax);
if (!converged) {
    std::cerr << "MGMRES did not converge.\n";
}
```

替换为：

```cpp
bool converged = false;

switch (matrixSolverType_) {
case 0:
    std::cout << "Matrix solver: GMRES (FMM)\n";
    converged = MGMRES::mgmres(
        n, x, rhs, itr_max, mr, tol_abs, tol_rel, ax);
    if (!converged) {
        std::cerr << "GMRES did not converge (FMM).\n";
    }
    break;

case 1:
    std::cout << "Matrix solver: CGS (FMM)\n";
    converged = MCGS::CGS::solve(
        n, x, rhs, itr_max, tol_rel, ax);
    if (!converged) {
        std::cerr
            << "CGS did not converge or encountered breakdown (FMM).\n";
    }
    break;

default:
    throw std::invalid_argument(
        "Invalid matrixSolverType_ in FMM "
        "(must be 0 for GMRES or 1 for CGS)");
}
```

不要把 `ax` 移到 `switch` 的某个分支内。两种迭代算法必须调用同一个 FMM 矩阵－向量乘过程，才能保证它们求解的是同一个线性系统。

---

## 4. 修改 `MLFMM.h`

### 4.1 构造函数增加参数

将原声明：

```cpp
MLFMM(
    const RCSExportConfig& cfg, const int selectIntegralEqu,
    const std::string selectMono_Dual, const std::string pol_wave,
    const std::vector<std::vector<OCTree::Node*>>& octreeNodes,
    const std::vector<std::vector<OCTree::nodePoint>>& octreeNodesDRvec_,
    const std::vector<RWGBase>& rwgs, const int maxLevel_,
    const gaussPoints& gausspoint, const double E0,
    const EMSource& wave);
```

改为：

```cpp
MLFMM(
    const RCSExportConfig& cfg, const int selectIntegralEqu,
    const int selectMatrixSolver,
    const std::string selectMono_Dual, const std::string pol_wave,
    const std::vector<std::vector<OCTree::Node*>>& octreeNodes,
    const std::vector<std::vector<OCTree::nodePoint>>& octreeNodesDRvec_,
    const std::vector<RWGBase>& rwgs, const int maxLevel_,
    const gaussPoints& gausspoint, const double E0,
    const EMSource& wave);
```

### 4.2 保存选择值

在 `integralEquType_` 后加入：

```cpp
const int integralEquType_;
const int matrixSolverType_; // 0 = GMRES, 1 = CGS
```

---

## 5. 修改 `MLFMM.cpp`

### 5.1 引入 CGS

在头文件区加入：

```cpp
#include <stdexcept>

#include "CGS.h"
```

即让相关部分包含：

```cpp
#include "GMRES.h"
#include "CGS.h"
```

### 5.2 修改构造函数定义

将构造函数头部和初始化列表改为：

```cpp
MLFMM::MLFMM(
    const RCSExportConfig& cfg,
    const int selectIntegralEqu,
    const int selectMatrixSolver,
    const std::string selectMono_Dual,
    const std::string pol_wave,
    const std::vector<std::vector<OCTree::Node*>>& octreeNodes,
    const std::vector<std::vector<OCTree::nodePoint>>& octreeNodesDRvec_,
    const std::vector<RWGBase>& rwgs,
    const int maxLevel_,
    const gaussPoints& gausspoint,
    const double E0,
    const EMSource& wave)
    : octreeNodes_(octreeNodes),
      octreeNodesDRvec_(octreeNodesDRvec_),
      rwgs(rwgs),
      maxLevel_(maxLevel_),
      E0(E0),
      integralEquType_(selectIntegralEqu),
      matrixSolverType_(selectMatrixSolver),
      gausspoint(gausspoint),
      wave(wave)
{
    // 原有构造函数主体保持不变
```

这里按 `MLFMM.h` 中的成员声明顺序书写初始化列表，可减少 MSVC 的初始化顺序警告。未显式初始化的 `L_k1`、`L_k2`、`levelSpan`、`row` 仍由后续 `computeBase()` 写入。

### 5.3 替换 `matrix_solver()` 末尾的固定调用

保留现有 `ax` lambda 的全部 MLFMA 上行、转移、下行和近场代码。仅把函数末尾：

```cpp
bool converged = MGMRES::mgmres(
    n, x, rhs, itr_max, mr, tol_abs, tol_rel, ax);
if (!converged) {
    std::cerr << "MGMRES did not converge.\n";
}
```

替换为：

```cpp
bool converged = false;

switch (matrixSolverType_) {
case 0:
    std::cout << "Matrix solver: GMRES (MLFMA)\n";
    converged = MGMRES::mgmres(
        n, x, rhs, itr_max, mr, tol_abs, tol_rel, ax);
    if (!converged) {
        std::cerr << "GMRES did not converge (MLFMA).\n";
    }
    break;

case 1:
    std::cout << "Matrix solver: CGS (MLFMA)\n";
    converged = MCGS::CGS::solve(
        n, x, rhs, itr_max, tol_rel, ax);
    if (!converged) {
        std::cerr
            << "CGS did not converge or encountered breakdown (MLFMA).\n";
    }
    break;

default:
    throw std::invalid_argument(
        "Invalid matrixSolverType_ in MLFMM "
        "(must be 0 for GMRES or 1 for CGS)");
}
```

---

## 6. 同步修改 `main.cpp` 中 FMM/MLFMM 的构造调用

你目前只把 `selectMatrixSolver` 传给了 `MoM`。完成上述构造函数修改后，还必须同步修改另外两处调用。

### 6.1 FMM

原代码：

```cpp
FMM fmm(cfg, selectIntegralEqu, selectMono_Dual, polarization,
    octree.octreeNodes_, rwgBases, octree.maxLevel_, gausspoint,
    E0, wave);
```

改为：

```cpp
FMM fmm(cfg, selectIntegralEqu, selectMatrixSolver,
    selectMono_Dual, polarization,
    octree.octreeNodes_, rwgBases, octree.maxLevel_, gausspoint,
    E0, wave);
```

### 6.2 MLFMM

原代码：

```cpp
MLFMM mlfmm(cfg, selectIntegralEqu, selectMono_Dual, polarization,
    octree.octreeNodes_, octree.octreeNodesDRvec_, rwgBases,
    octree.maxLevel_, gausspoint, E0, wave);
```

改为：

```cpp
MLFMM mlfmm(cfg, selectIntegralEqu, selectMatrixSolver,
    selectMono_Dual, polarization,
    octree.octreeNodes_, octree.octreeNodesDRvec_, rwgBases,
    octree.maxLevel_, gausspoint, E0, wave);
```

你现在对 `selectMatrixSolver` 的合法性校验位于三种算法分支之前，因此 MoM、FMM、MLFMM 都能共享该校验，无需在每个分支重复。

---

## 7. 修改 `BatchRunner.h`

BatchRunner 需要完成“配置文件 → 批配置 → 单个任务 → 求解对象”的完整传递。少任意一层，配置值都会丢失。

### 7.1 `SimTask` 增加字段

在 `selectIntegralEqu` 后加入：

```cpp
int selectAlgorithm;       // 0 = MoM, 1 = FMM, 2 = MLFMM
int selectIntegralEqu;     // 0 = EFIE, 1 = CFIE, 2 = PMCHWT
int selectMatrixSolver;    // 0 = GMRES, 1 = CGS
std::string selectMono_Dual;
```

### 7.2 `BatchConfig` 增加默认值

在相应配置区域加入：

```cpp
int selectAlgorithm = 2;       // 默认 MLFMM
int selectIntegralEqu = 1;     // 默认 CFIE
int selectMatrixSolver = 0;    // 默认 GMRES
std::string selectMono_Dual = "dual";
```

默认使用 GMRES 可以保持旧批处理配置文件的行为：旧文件即使没有 `matrix_solver` 项，也仍然选择 GMRES。

### 7.3 `generateTasks()` 复制字段

在生成任务时加入：

```cpp
task.selectAlgorithm = selectAlgorithm;
task.selectIntegralEqu = selectIntegralEqu;
task.selectMatrixSolver = selectMatrixSolver;
task.selectMono_Dual = selectMono_Dual;
```

如果漏掉这一行，解析器虽然读到了配置值，但每个 `SimTask` 中的值会处于未初始化状态。

### 7.4 解析 `matrix_solver`

在 `BatchConfigParser::parse()` 中，紧跟 `integral_equation` 分支加入：

```cpp
else if (key == "integral_equation") {
    cfg.selectIntegralEqu = std::stoi(value);
}
else if (key == "matrix_solver") {
    cfg.selectMatrixSolver = std::stoi(value);
}
else if (key == "mono_dual") {
    cfg.selectMono_Dual = value;
}
```

### 7.5 在任务日志中显示求解器

把：

```cpp
std::cout << "  Algo: " << task.selectAlgorithm
    << ", IE: " << task.selectIntegralEqu << "\n";
```

改为：

```cpp
std::cout << "  Algo: " << task.selectAlgorithm
    << ", IE: " << task.selectIntegralEqu
    << ", Solver: " << task.selectMatrixSolver
    << " (0=GMRES, 1=CGS)\n";
```

### 7.6 尽早校验选择值

在 `executeSingleTask()` 的 `try` 块开始处、读取网格之前加入：

```cpp
if (task.selectMatrixSolver != 0 && task.selectMatrixSolver != 1) {
    throw std::invalid_argument(
        "Invalid matrix_solver value: "
        + std::to_string(task.selectMatrixSolver)
        + " (must be 0 for GMRES or 1 for CGS)");
}
```

放在读取网格和构建八叉树之前，可以避免无效配置消耗大量预处理时间。

这里先输出数值而不是用“非 0 即 CGS”的三元表达式。任务摘要发生在 `executeSingleTask()` 校验之前；直接输出数值能避免非法值（例如 `2`）被错误显示成 CGS。

### 7.7 修改三种求解对象的构造调用

把 `BatchRunner.h` 中的整个“执行求解”分支改为：

```cpp
// 执行求解
if (task.selectAlgorithm == 0) {
    MoM mom(cfg, task.selectIntegralEqu, task.selectMatrixSolver,
        task.selectMono_Dual, pol,
        octree.octreeNodes_, rwgBases,
        octree.maxLevel_, gausspoint, E0, wave);
}
else if (task.selectAlgorithm == 1) {
    FMM fmm(cfg, task.selectIntegralEqu, task.selectMatrixSolver,
        task.selectMono_Dual, pol,
        octree.octreeNodes_, rwgBases,
        octree.maxLevel_, gausspoint, E0, wave);
}
else if (task.selectAlgorithm == 2) {
    MLFMM mlfmm(cfg, task.selectIntegralEqu, task.selectMatrixSolver,
        task.selectMono_Dual, pol,
        octree.octreeNodes_, octree.octreeNodesDRvec_,
        rwgBases, octree.maxLevel_, gausspoint, E0, wave);
}
```

其中 MoM 调用也必须增加 `task.selectMatrixSolver`。你已经修改了 MoM 构造函数，而当前 BatchRunner 仍使用旧参数列表；若不修改，这里会出现“没有匹配的重载函数”编译错误。

---

## 8. 修改 `batch_config.txt`

在仿真参数区域加入：

```ini
# 0 = GMRES, 1 = CGS
matrix_solver = 0
```

完整示例：

```ini
# ============ 仿真参数 ============
# algorithm: 0 = MoM, 1 = FMM, 2 = MLFMM
algorithm = 2

# integral_equation: 0 = EFIE, 1 = CFIE, 2 = PMCHWT
integral_equation = 1

# matrix_solver: 0 = GMRES, 1 = CGS
matrix_solver = 1

mono_dual = dual
gauss_points = 7
omp_threads = 14
```

此处 `matrix_solver = 1` 表示该批次生成的全部任务和两个极化都使用 CGS。若将其设为 `0`，则全部使用 GMRES。

---

## 9. 为什么不修改 `RCS.h`

你已经把 `RCS.h` 中四个求解位置统一改成：

```cpp
solver.matrix_solver(...);
```

只要 `MoM`、`FMM`、`MLFMM` 都提供同签名的 `matrix_solver()`，模板就会在编译时调用对应类的方法。算法选择由各对象内部保存的 `matrixSolverType_` 决定，所以 `RCS.h` 不需要知道 GMRES 或 CGS 的具体类型。

这种分层方式的职责比较清楚：

```text
RCS.h               负责“何时求解”
MoM/FMM/MLFMM       负责“A·x 如何计算”
matrixSolverType_   负责“使用哪种迭代过程”
GMRES.h / CGS.h     负责具体迭代算法
```

---

## 10. 参数语义

两种求解器共用 `matrix_solver()` 的形参，但并非所有参数都适用于 CGS：

| 参数 | GMRES | CGS |
|---|---|---|
| `itr_max` | 最大总迭代次数 | 最大迭代次数 |
| `mr` | Krylov 子空间/重启参数 | 不使用 |
| `tol_abs` | 绝对残差阈值 | 当前实现不使用 |
| `tol_rel` | 相对残差阈值 | 作为 `tol` 使用 |

当前 `CGS.h` 使用：

```text
RSS = ||b - A·x|| / ||b||
```

并以 `RSS < tol` 判断收敛，因此把 `tol_rel` 传给 CGS 是正确对应。不要把 `mr` 或 `tol_abs` 硬套进 CGS 调用。

---

## 11. FMM/MLFMA 下使用 CGS 的计算特点

一次 GMRES 内迭代通常调用一次 `ax`，并需要保存 Krylov 基；一次 CGS 迭代通常调用两次 `ax`，但不保存整个 Krylov 基。因此：

- CGS 一般占用更少的迭代器工作内存；
- 对 FMM/MLFMA 而言，`ax` 本身很昂贵，CGS 单次迭代的时间不一定更短；
- CGS 残差可能明显振荡，也可能发生 `rho` 或 `alpha` 分母接近零的 breakdown；
- GMRES 通常更稳健，但 `mr` 较大时需要更多内存。

所以应该根据总矩阵－向量乘次数、总耗时、最终残差和 RCS 误差共同比较，而不是只比较输出的迭代编号。

---

## 12. 推荐验证顺序

### 12.1 编译验证

先确认以下接口完全一致：

```text
FMM.h       ↔ FMM.cpp
MLFMM.h     ↔ MLFMM.cpp
main.cpp    ↔ 三个类的构造函数
BatchRunner ↔ 三个类的构造函数
RCS.h       ↔ 三个类的 matrix_solver()
```

### 12.2 单次运行验证

依次测试：

```cpp
selectAlgorithm = 1;
selectMatrixSolver = 0; // FMM + GMRES
```

```cpp
selectAlgorithm = 1;
selectMatrixSolver = 1; // FMM + CGS
```

```cpp
selectAlgorithm = 2;
selectMatrixSolver = 0; // MLFMA + GMRES
```

```cpp
selectAlgorithm = 2;
selectMatrixSolver = 1; // MLFMA + CGS
```

### 12.3 批处理验证

分别用：

```ini
matrix_solver = 0
```

和：

```ini
matrix_solver = 1
```

运行一个规模较小的任务，确认任务摘要与每次求解输出的名称一致。

### 12.4 数值验证

对同一模型、同一积分方程和同一误差阈值，比较：

1. 最终相对残差；
2. GMRES 与 CGS 的矩阵－向量乘次数；
3. 总求解时间；
4. 电流系数或最终 RCS；
5. EFIE、CFIE、PMCHWT 三种系统是否都能运行。

GMRES 与 CGS 的解不必逐位完全相同，但在都正常收敛时，最终物理量的差异应与迭代误差量级相符。

---

## 13. 必要修改文件汇总

```text
FMM.h              构造参数 + matrixSolverType_ 成员
FMM.cpp            包含 CGS.h + 初始化成员 + 求解器分派
MLFMM.h            构造参数 + matrixSolverType_ 成员
MLFMM.cpp          包含 CGS.h + 初始化成员 + 求解器分派
main.cpp           向 FMM/MLFMM 传递 selectMatrixSolver
BatchRunner.h      配置、任务、解析、校验和三种构造调用
batch_config.txt   增加 matrix_solver 配置项
```

无需继续修改：

```text
RCS.h
GMRES.h
CGS.h
MoM.h
MoM.cpp
```

本文档只提供具体修改代码和解释，没有修改任何现有 C++ 源文件或批处理配置文件。
