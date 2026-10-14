# HSB-MLFMA-Moving 代码修改方案

## 1. 算法定位

`HSB-MLFMA-Moving` 定义为：

```text
运动目标 + 单站 RCS + MLFMA + HSB 路径角域压缩
```

第一版只考虑单站雷达，不考虑双站。发射和接收始终为同一位置：

```text
Tx = Rx = radar_position
```

目标是通过运动路径反推出每个时刻在目标本体系中的单站角度，只为这些角度附近生成 HSB 准确采样点，从而减少运动目标 RCS 时间序列的计算量。

## 2. 第一版支持范围

第一版建议限定为：

```text
目标类型：PEC
算法：MLFMA
收发模式：Mono
积分方程：优先 CFIE，EFIE 可兼容
姿态：刚体旋转
位置：刚体平移
路径输入：txt 文件
RCS 输出：运动时间序列 RCS
```

物理假设：

1. 目标几何不变，RWG 基函数不变。
2. 目标为 PEC 刚体，阻抗矩阵 `Z` 不随目标平移或旋转改变。
3. 目标旋转等效为入射/观测方向反向旋转到目标本体系。
4. 单站 RCS 中发射方向和接收方向绑定，不引入独立双站角。
5. 平移不改变 RCS 幅值，只改变雷达到目标的视线方向和距离。

第一版不建议加入距离相位。后续若需要复回波或微多普勒，再单独加入：

```text
E_echo(t) = E_body(t) * exp(-j * 2 * k * range(t))
```

## 3. 运动路径 txt 格式

建议第一版路径文件只支持固定四元数格式：

```text
# time_s  x  y  z  qw  qx  qy  qz
0.000     1000 0   0   1.0       0.0 0.0 0.0
0.010     1000 5   0   0.9999619 0.0 0.0 0.0087265
0.020     1000 10  0   0.9998477 0.0 0.0 0.0174524
```

列含义：

```text
time_s        时间，单位 s
x y z         目标质心在全局坐标系的位置
qw qx qy qz   目标姿态四元数，表示 body -> global 的旋转
```

单站雷达位置不写入每一行，而是放在配置中：

```text
radar_position = x,y,z
```

## 4. 角度映射关系

每个时刻从路径文件读取：

```text
p_target(t)
q_body_to_global(t)
```

单站雷达位置：

```text
p_radar
```

目标指向雷达的全局观察方向：

```text
r_global(t) = normalize(p_radar - p_target(t))
```

四元数转旋转矩阵：

```text
R(t): body -> global
```

映射到目标本体系：

```text
r_body(t) = R(t)^T * r_global(t)
```

再转成 HSB 查询角：

```text
theta_body = acos(r_body.z)
phi_body   = atan2(r_body.y, r_body.x)
```

角度范围：

```text
theta_body_deg in [0, 180]
phi_body_deg   in [0, 360)
```

极点规则沿用当前 RCS 输出优化方案：

```text
theta = 0 或 180 时，phi 统一视为 0
```

## 5. RCSExportConfig 修改

建议在 `Code/RCSExportConfig.h` 的 `RCSExportConfig` 中新增：

```cpp
bool useHSBMLFMAMoving{ false };
std::string movingPathFile;

double radarX{ 0.0 };
double radarY{ 0.0 };
double radarZ{ 0.0 };

double movingAngleMarginDeg{ 2.0 };
double movingAngleMergeTolDeg{ 0.25 };

bool exportMovingAngles{ true };
bool exportMovingComplexField{ false };
bool exportMovingSamples{ false };
```

字段含义：

```text
useHSBMLFMAMoving        启用 HSB-MLFMA-Moving
movingPathFile           运动路径 txt 文件
radarX/Y/Z               单站雷达全局坐标，Tx=Rx
movingAngleMarginDeg     路径角域压缩时的安全裕度
movingAngleMergeTolDeg   查询角去重阈值
exportMovingAngles       输出路径角度映射检查文件
exportMovingComplexField 输出 Ev/Eh 复数
exportMovingSamples      输出实际 HSB 准确采样点
```

`validate()` 中建议增加：

```text
if useHSBMLFMAMoving:
    selectAlgo 必须为 MLFMA
    material 必须为 PEC
    monoDual 必须为 Mono
    integralEqu 不应为 PMCHWT
    movingPathFile 非空
    movingAngleMarginDeg >= 0
    movingAngleMergeTolDeg > 0
```

## 6. 批处理配置修改

在 `Code/BatchRunner.h` 的 `BatchConfig` 和 `SimTask` 中增加同名字段，并在 parser 增加 key。

配置文件示例：

```text
select_algorithm = 2
material = pec
integral_equation = 1
mono_dual = mono
use_hsb_mlfma_moving = true

moving_path_file = D:\MyCode\PMCHWT_MLFMA\DATA\moving_path.txt
radar_position = 0,0,0

hsb_rho = 1.2
moving_angle_margin_deg = 2.0
moving_angle_merge_tol_deg = 0.25

export_moving_angles = true
export_moving_complex_field = false
export_moving_samples = false
```

保留现有 `use_hsb_mono`，但建议优先级为：

```text
use_hsb_mlfma_moving > use_hsb_mono > 普通 mono
```

当 `use_hsb_mlfma_moving=true` 时，不再走普通全角域 HSB 输出。

## 7. 新增 MovingTarget 模块

建议新增：

```text
Code/MovingTarget.h
Code/MovingTarget.cpp
```

用于运动路径读取、四元数处理、坐标变换和角度去重。

建议数据结构：

```cpp
struct MovingState {
    double time_s = 0.0;
    double position[3] = { 0.0, 0.0, 0.0 };
    double quat[4] = { 1.0, 0.0, 0.0, 0.0 };
};

struct MovingAngleSample {
    double time_s = 0.0;

    double pos[3] = { 0.0, 0.0, 0.0 };
    double range = 0.0;

    double thetaGlobalDeg = 0.0;
    double phiGlobalDeg = 0.0;

    double thetaBodyDeg = 0.0;
    double phiBodyDeg = 0.0;
};
```

建议函数命名：

```cpp
namespace MovingTarget {
    std::vector<MovingState> readMovingPathFile(const std::string& path);

    void normalizeQuaternion(double q[4]);

    void quaternionToMatrix(const double q[4], double R[3][3]);

    std::array<double, 3> normalizeVector(const std::array<double, 3>& v);

    std::array<double, 3> rotateGlobalToBody(
        const double R[3][3],
        const std::array<double, 3>& vGlobal);

    void vectorToThetaPhiDeg(
        const std::array<double, 3>& v,
        double& thetaDeg,
        double& phiDeg);

    MovingAngleSample computeMovingAngleSample(
        const MovingState& state,
        const std::array<double, 3>& radarPos);

    std::vector<MovingAngleSample> computeMovingAngleSamples(
        const std::vector<MovingState>& states,
        const std::array<double, 3>& radarPos);

    std::vector<MovingAngleSample> mergeMovingAngleSamples(
        const std::vector<MovingAngleSample>& samples,
        double mergeTolDeg);
}
```

路径文件读取规则：

1. 空行跳过。
2. `#` 开头的注释行跳过。
3. 每行必须包含 8 列：`time_s x y z qw qx qy qz`。
4. 四元数自动归一化。
5. 如果目标位置与雷达位置距离过小，抛出异常。
6. `theta/phi` 转换后归一化到合法范围。

## 8. 路径角域压缩

输入为所有时刻的目标本体系查询角：

```text
movingAngles = { thetaBodyDeg(t_i), phiBodyDeg(t_i) }
```

先做查询角去重：

```text
uniqueMovingAngles = mergeMovingAngleSamples(movingAngles, movingAngleMergeTolDeg)
```

去重判断：

```text
abs(theta1 - theta2) <= tol
deltaPhi(phi1, phi2) <= tol
```

其中：

```cpp
deltaPhi = min(abs(phi1 - phi2), 360 - abs(phi1 - phi2))
```

极点附近：

```text
theta = 0 或 180 时，phi 统一为 0
```

## 9. HSB 准确采样点生成

当前 `computeMonoStatic_HSB` 的全局建库方式为：

```text
theta_m  = pi * m / (thetaCount + 1)
phiOrder = ceil(W * sin(theta_m))
phi_n    = 360 * n / (2 * phiOrder + 1)
```

`HSB-MLFMA-Moving` 不改变这套 HSB 采样公式，而是做路径驱动筛选：

1. 根据 `W` 得到理论 `thetaCount` 和每个 ring 的 `phiOrder`。
2. 对每个 body-frame 查询角 `theta_q / phi_q`，选择附近的若干 `theta ring`。
3. 对每个被选中的 `theta ring`，选择 `phi_q` 附近的若干 `phi_n`。
4. 汇总并去重，形成稀疏 HSB `rings`。

建议新增在 `Code/OverloadAlgo.h` 的 `namespace my_hsb`：

```cpp
inline std::vector<int> selectMovingThetaRingIndices(
    double thetaQueryDeg,
    int thetaCount,
    double marginDeg);

inline std::vector<int> selectMovingPhiIndices(
    double phiQueryDeg,
    int phiOrder,
    double marginDeg);

inline std::vector<std::vector<HSBMonoSample>> buildMovingHSBSampleLayout(
    const std::vector<MovingAngleSample>& queryAngles,
    double W,
    double marginDeg);
```

`buildMovingHSBSampleLayout` 只生成角度布局：

```text
sample.thetaDeg
sample.phiDeg
sample.phiOrder
sample.field 为空
```

真正求解 `field` 仍放在 RCS 求解流程中。

## 10. 不使用矩形角域的原因

如果运动路径穿过：

```text
phi = 350 -> 10
```

矩形角域容易错误变成：

```text
phi = 10..350
```

这会几乎覆盖整个方位角，计算量暴涨。

路径点集筛选方式是：

```text
围绕每个实际查询角选邻域采样点
统一去重
```

它能自然处理 `0/360` 周期边界，更适合运动路径。

## 11. RCS.h 新增主函数

建议在 `Code/RCS.h` 中新增：

```cpp
template<typename SolverType, typename ComputeVFunc>
inline void computeMonoStatic_HSB_Moving(
    SolverType& solver,
    const RCSExportConfig& cfg,
    ComputeVFunc computeV,
    double rho = 1.2)
```

函数流程：

```text
1. 读取 movingPathFile。
2. 计算每个时刻的 MovingAngleSample。
3. 若 exportMovingAngles=true，输出角度映射检查文件。
4. 对 body-frame 查询角去重。
5. 根据目标包围球和波数计算 W。
6. 用 uniqueMovingAngles 生成稀疏 HSB rings。
7. 对 rings 中每个 HSB 采样点执行准确求解：
   a. pw.initPW(data_pol, thetaDeg, phiDeg)
   b. computeV(solver, kinc, einc, hinc)
   c. solver.matrix_solver(...)
   d. calculateRCS(...)
   e. sample.field = RCS_pre
8. 若 exportMovingSamples=true，输出实际 HSB 准确采样点。
9. 对原始每个运动时刻：
   a. 用 thetaBodyDeg / phiBodyDeg 调用 interpolateHSBField
   b. 用 projectRCSGridValue 投影出 RCS
   c. 写运动 RCS 时间序列
```

注意：最终输出使用原始全部时刻，不使用去重后的时刻。去重只用于减少 HSB 准确求解点。

## 12. 输出文件设计

运动 RCS 输出：

```text
*_RCS_HSB_MLFMA_MOVING_*.txt
```

列：

```text
idx
time_s
range_m
theta_global_deg
phi_global_deg
theta_body_deg
phi_body_deg
rcs_db
rcs_v_db
rcs_h_db
```

角度映射检查输出：

```text
*_ANGLES_HSB_MLFMA_MOVING_*.txt
```

列：

```text
idx
time_s
pos_x
pos_y
pos_z
range_m
theta_global_deg
phi_global_deg
theta_body_deg
phi_body_deg
```

可选复场输出：

```text
*_FIELD_HSB_MLFMA_MOVING_*.txt
```

列：

```text
idx
time_s
range_m
theta_body_deg
phi_body_deg
Ev_re
Ev_im
Eh_re
Eh_im
```

可选 HSB 样本输出：

```text
*_FIELD_HSB_MLFMA_MOVING_SAMPLES_*.txt
```

列格式可复用现有 `exportHSBSamples` 的格式。

## 13. MLFMM.cpp 接入

建议在 `Code/MLFMM.cpp` 的 PEC 单站分支接入。

重点接入位置：

```cpp
void MLFMM::mlfmm_Mono_Pec_Cfie(
    const RCSExportConfig& cfg,
    const std::string pol_wave)
```

建议优先级：

```cpp
if (cfg.useHSBMLFMAMoving) {
    RCSUtils::computeMonoStatic_HSB_Moving(*this, cfg, computeV, cfg.hsbRho);
}
else if (cfg.useHSBMono) {
    RCSUtils::computeMonoStatic_HSB(*this, cfg, computeV, cfg.hsbRho);
}
else {
    RCSUtils::computeMonoStatic(*this, cfg, computeV);
}
```

`mlfmm_Mono_Pec_Efie` 也可以同样接入，但第一版建议先接 CFIE，确认无误后再扩展 EFIE。

## 14. main.cpp 单次运行参数

建议在 `Code/main.cpp` 中增加本地调试变量：

```cpp
bool useHSBMLFMAMoving = false;
std::string movingPathFile = "D:\\MyCode\\PMCHWT_MLFMA\\DATA\\moving_path.txt";

double radarX = 0.0;
double radarY = 0.0;
double radarZ = 0.0;

double movingAngleMarginDeg = 2.0;
double movingAngleMergeTolDeg = 0.25;

bool exportMovingAngles = true;
bool exportMovingComplexField = false;
bool exportMovingSamples = false;
```

并赋值给 `cfg`：

```cpp
cfg.useHSBMLFMAMoving = useHSBMLFMAMoving;
cfg.movingPathFile = movingPathFile;
cfg.radarX = radarX;
cfg.radarY = radarY;
cfg.radarZ = radarZ;
cfg.movingAngleMarginDeg = movingAngleMarginDeg;
cfg.movingAngleMergeTolDeg = movingAngleMergeTolDeg;
cfg.exportMovingAngles = exportMovingAngles;
cfg.exportMovingComplexField = exportMovingComplexField;
cfg.exportMovingSamples = exportMovingSamples;
```

调试时推荐：

```text
selectAlgorithm = 2
selectMono_Dual = "mono"
selectIntegralEqu = 1
useHSBMLFMAMoving = true
useHSBMono = false
```

## 15. Visual Studio 项目文件修改

因为工程使用 Visual Studio，需要把新增文件加入：

```text
Code/MLFMA.vcxproj
Code/MLFMA.vcxproj.filters
```

新增编译单元：

```text
MovingTarget.cpp
```

新增头文件：

```text
MovingTarget.h
```

## 16. 代码复用点

应尽量复用现有函数：

```text
RCSUtils::calculateRCS
RCSUtils::projectRCSGridValue
RCSUtils::isThetaPole
RCSUtils::canonicalPolePhiDeg
my_hsb::estimateBoundingSphereRadius
my_hsb::normalizePhiDeg
my_hsb::interpolateHSBField
RCSUtils::exportHSBSamples
solver.matrix_solver
MLFMM 分支中的 computeV lambda
```

新增代码只负责：

```text
运动路径读取
姿态/坐标变换
路径角度去重
HSB 采样点筛选
运动时间序列输出
```

## 17. 验证方案

### 17.1 路径角度映射验证

目标不旋转：

```text
q = 1,0,0,0
```

目标沿直线移动。检查：

```text
theta_global == theta_body
phi_global == phi_body
```

### 17.2 固定目标姿态验证

所有时刻位置相同、姿态相同。输出 RCS 应完全相同。

### 17.3 纯旋转验证

目标位置固定，姿态绕 z 轴旋转。预期：

```text
theta_body 基本不变
phi_body 随姿态反向变化
```

### 17.4 与直接单站计算对比

抽取 5-10 个时刻，用普通 MLFMA mono 直接计算对应 `theta_body/phi_body`，与 `HSB-MLFMA-Moving` 插值结果对比。

重点调参：

```text
hsb_rho
movingAngleMarginDeg
movingAngleMergeTolDeg
```

## 18. 建议实现顺序

1. 增加 `RCSExportConfig` 配置字段。
2. 增加 `MovingTarget.h/.cpp`，实现路径读取和角度映射。
3. 增加角度映射输出文件，先不接 HSB。
4. 实现路径角度去重。
5. 在 `my_hsb` 中增加路径驱动 HSB 样本布局生成函数。
6. 新增 `computeMonoStatic_HSB_Moving`。
7. 接入 `mlfmm_Mono_Pec_Cfie`。
8. 增加 batch 配置解析。
9. 加入 `MLFMA.vcxproj` 和 `MLFMA.vcxproj.filters`。
10. 编译验证。
11. 用小路径文件做数值验证。

## 19. 需要审核确认的设计选择

建议第一版固定以下规则：

```text
1. HSB-MLFMA-Moving 只支持单站，不支持双站。
2. 第一版只支持 PEC，不支持 PMCHWT 介质目标。
3. 运动路径 txt 第一版只支持 time x y z qw qx qy qz。
4. 雷达为单站点：Tx=Rx=radar_position。
5. 第一版只输出 RCS 幅值，不乘平移相位。
6. HSB 样本使用路径点集邻域筛选，而不是 theta/phi 矩形角域。
7. 第一版优先接入 MLFMA + PEC + Mono + CFIE。
```

这套设计的主线是：先通过运动路径得到目标本体系中的单站角度，再只为这些角度附近建立 HSB 样本库。这样可以最大化复用现有 MLFMA/HSB 代码，同时实质性减少运动目标时间序列 RCS 的计算量。
