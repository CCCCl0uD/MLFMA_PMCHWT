# GMRES / CGS OpenMP 验证入口

在项目根目录的 PowerShell 中运行：

```powershell
.\tests\run_solver_openmp.ps1
.\tests\run_solver_openmp.ps1 -Benchmark -Size 100000 -Threads 4
```

需要已安装 MSVC C++ 工具和 Windows SDK。脚本通过 `vswhere` 找到编译器，将构建和日志写入已忽略的 `output/solver_openmp/`，不会修改主工程的线程配置或现有求解器文件。

默认参考版本固定为改造前提交 `9a357cc7304f85439bf1d85926e3b9e9f866cb91`，通过 Git 原始字节流导出其两个头文件，并重命名命名空间，使新版与原版可在同一进程中比较。若仓库为浅克隆，须先确保该提交存在。也可通过 `-BaselineRef <commit>` 指定兼容的原始参考；参考头文件应为改造前、可独立包含的版本。

验证程序覆盖：

- 共轭/非共轭内积、范数、残差与影子残差、融合残差更新、分块解更新。
- 阈值前后、不整除线程数的长度、小向量回退和外层并行区内的单调用者回退。
- 1/2/4/8/12/16 线程，每档重复 3 次；同时检查长向量实际使用了多个工作线程。
- 已知解的复数非厄米对角占优系统，零/非零初值，GMRES 重启、CGS 无预处理与对角预处理。
- 原始版本的返回状态、解向量与 `ax` 调用次数；独立计算真实残差，而非只相信递推残差。
- 零右端、初始收敛、迭代次数耗尽、CGS 击穿及 GMRES 原有线性相关处理。
- OpenMP、仅关闭求解器并行、完全不启用 OpenMP 三种编译方式。

无 OpenMP 的验证仅针对两个求解器头文件，不表示整个主工程已支持无 OpenMP 构建。当前回归覆盖成功算例 468 组（另有对应原始版本求解及边界/内核检查）。测试失败时返回非零退出码。

基准程序参数顺序：

```powershell
.\output\solver_openmp\benchmark_openmp.exe 100000 4 5 50
# n、线程数、测量次数、GMRES 重启维度
```

CSV 中 `current` 为新版，`baseline` 为原始串行求解器；两者的矩阵向量乘均使用同样的线程数。`GMRES/CGS` 行给出完整合成求解时间、累计 `ax` 时间的独立中位数、调用次数与真实残差；`dot/MGS2` 行给出内积和两轮修正 Gram–Schmidt 的微基准。计时前预热，计时区外重新计算真实残差。`ax_seconds` 与总时间分别取中位数，不保证来自同一次样本。

每个基准的矩阵是 O(n) 的三对角复数矩阵；它用于隔离求解器内部向量工作的收益，不能代替 MoM/FMM/MLFMM 的端到端计时。生成基向量的时间不计入 MGS 微基准。运行时应避免同时进行其他编译或仿真，以降低测量噪声。

编译控制项：

| 宏 | 默认值/用途 |
| --- | --- |
| `MLFMA_SOLVER_PARALLEL_THRESHOLD` | `32768`，按完整未知量向量的复数元素数量判定；PMCHWT 为 `2*row` |
| `MLFMA_SOLVER_DISABLE_OPENMP` | 默认未定义；定义后仅关闭求解器向量内核并行，保留 `ax` 原有并行 |

这些宏应对包含求解器的所有翻译单元保持一致，避免头文件内联定义不一致。更改后重新编译，不能仅通过设置环境变量改变阈值。线程数仍由 `main.cpp` 或批处理 `omp_threads` 管理。

真实工程算例的 RCS 文件可单独配对检查：

```powershell
.\tests\compare_solver_rcs.ps1 `
  -ParallelDirectory output\solver_openmp\engineering\parallel `
  -BaselineDirectory output\solver_openmp\engineering\baseline
```

该脚本递归匹配对应配置目录中除时间戳外同名的 RCS 文件，要求每个文件恰有一个参考版本。默认要求各极化通道线性散射功率差除以参考总功率后不超过 `1e-4`，同时报告总 RCS 和任意通道的最大 dB 差，避免只用弱交叉极化分量的 dB 差判定整体误差。也可通过 `-RelativePowerTolerance` 明确指定所需精度。该检查不代替真实残差或电流解验证。
