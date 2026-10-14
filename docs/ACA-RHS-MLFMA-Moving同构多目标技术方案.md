# ACA-RHS-MLFMA-Moving 同构多目标技术方案

## 1. 方案结论

本方案把一份 NAS 描述的单目标视为几何模板，由程序按照多个初始位置生成同构目标实例，再把所有实例组成一个复合散射体：

```text
一份 NAS 模板 + 程序内的多个初始位置
                ↓
      同构多目标复合几何
                ↓
一个全局 OCTree + 一个固定 MLFMA 阻抗算子
                ↓
      ACA 压缩复合 RHS 时间序列
                ↓
逐时刻合成全体电流并计算相干总 RCS
```

“多目标阻抗矩阵”必须是包含完整互耦的联合系统：

```text
| Z00  Z01  ...  Z0,M-1 | | I0 |   | b0 |
| Z10  Z11  ...  Z1,M-1 | | I1 | = | b1 |
| ...                    | | .. |   | .. |
| ZM-1,0       ZM-1,M-1 | | IM |   | bM |
```

- `Zaa` 是实例 `a` 的自作用块；同构目标仅平移时，各对角块数学上相同。
- `Zab (a != b)` 是目标间互耦块，取决于相对位置，不能省略。
- MLFMA 不显式存储完整稠密矩阵：近区进入 `Z_near`，远区由聚合、转移和解聚隐式计算。
- 第一阶段先保证结果与“人工把多个目标合并成一个 NAS”一致，再加入同构复用优化。

推荐首版范围：

```text
算法：MLFMA
材料：PEC
方程：EFIE / CFIE
模式：Mono
网格：仅输入一份单目标 NAS
实例：相同网格、尺度、材料和朝向，仅初始位置不同
运动：所有实例保持相对位置不变，作为一个刚性编队整体运动
耦合：保留所有目标内自作用和目标间电磁互耦
加速：ACA-RHS + 同构几何/RHS/近场积分/远场积分复用
```

## 2. 物理边界

### 2.1 首版采用刚性编队运动

在编队坐标系中，实例 `a` 的初始位置为 `d_a`。运动文件的位置和四元数描述整个编队参考点：

```text
p_a_global(t) = p_formation(t) + R_formation(t) * d_a
R_a_global(t) = R_formation(t)
```

各实例相对位置和姿态始终不变，因此阻抗算子固定，现有形式仍成立：

```text
Z * I(t) = b(t)
```

复合几何、OCTree、`Z_near` 和 MLFMA 转移数据均只建立一次。

### 2.2 首版不支持各目标独立运动

若各目标拥有独立轨迹，则相对位移随时间变化，目标间互耦块也随时间变化：

```text
Zab = Zab(t)
```

这时不能继续使用“固定 `Z`、仅压缩 RHS”的 ACA-RHS-Moving。需要逐时刻更新算子或另行开发时变低秩互耦模型。因此首版必须明确校验：所有目标是一个刚性编队，而不是独立运动体。

### 2.3 远场照射条件

当前 Moving 实现用一条雷达视线生成目标本体系平面波方向。多目标编队继续使用这一模型，因此雷达对整个编队应满足远场条件。设编队最大尺寸为 `D_array`，程序应输出：

```text
Fraunhofer distance ≈ 2 * D_array^2 / lambda
```

若雷达距离不足，实例间的入射方向和距离相位不能用单一平面波准确描述，需要球面波模型；首版只警告，不扩展该能力。

## 3. 当前代码基础与接入点

| 现有能力 | 位置 | 多目标作用 |
|---|---|---|
| NAS 点和三角面读取 | `Code/NASReader.cpp` | 只读取一次模板 |
| RWG 生成 | `Code/RWGGenerator.cpp` | 生成模板 RWG 拓扑 |
| 全局 OCTree | `Code/OCTree.cpp` | 对复合几何建一棵树 |
| 近场阻抗 | `Code/Zmartix.h` | 计算目标内和目标间近区作用 |
| MLFMA 远区 | `Code/MLFMM.cpp`、`Code/ProcessMLFMM.h` | 隐式计算全部远区互耦 |
| RHS 单元素 | `Code/Vm.h` | ACA 行扫描与模板 RHS 复用 |
| ACA-RHS | `Code/ACARHS.h` | 压缩复合 RHS 矩阵 |
| 运动角度 | `Code/MovingTarget.cpp` | 编队整体运动映射到本体系角度 |
| RCS | `Code/RCS.h` | 合成所有实例的相干散射场 |

当前主链路是：

```text
readNasData → GenerateRWG → OCTree → MLFMM
```

多目标构造器插入在 `GenerateRWG` 与 `OCTree` 之间，把单目标模板扩展成连续的全局 `rwgs`，下游仍只面对一个联合系统。

## 4. 坐标与初始位置

NAS 坐标定义为模板局部坐标 `r_template`。增加模板参考点：

```cpp
std::array<double, 3> templateReferencePoint{ 0.0, 0.0, 0.0 };
```

实例坐标为：

```text
r_a,i = r_template_i - templateReferencePoint + initialPosition[a]
```

所以 `initialPosition[a]` 表示模板参考点在编队坐标系中的位置。首版建议 NAS 已以目标质心或安装基准点为原点。不能擅自用包围盒中心代替参考点，因为它未必是物理质心。

为避免后续相位公式混淆，定义实际施加到 NAS 原坐标上的平移量：

```text
delta_a = initialPosition[a] - templateReferencePoint
```

如果模板参考点为原点，则 `delta_a == initialPosition[a]`。

建议在 `main.cpp` 参数区给出位置：

```cpp
bool useIdenticalMultiTarget = true;

std::array<double, 3> templateReferencePoint = {
    0.0, 0.0, 0.0
};

std::vector<TargetInstanceConfig> targetInstances = {
    { 0, { 0.0,  0.0, 0.0 } },
    { 1, { 5.0,  0.0, 0.0 } },
    { 2, {-5.0,  0.0, 0.0 } }
};
```

首版只开放平移，不开放实例间不同旋转，以确保自阻抗、模板 RHS 和模板远场核可以精确复用。

## 5. 数据结构

建议新增：

```text
Code/IdenticalTargetAssembly.h
Code/IdenticalTargetAssembly.cpp
```

核心结构：

```cpp
struct TargetInstanceConfig {
    int instanceId = -1;
    std::array<double, 3> initialPosition{ 0.0, 0.0, 0.0 };
};

struct IdenticalTargetLayout {
    int targetCount = 0;
    int pointsPerTarget = 0;
    int facesPerTarget = 0;
    int rwgsPerTarget = 0;
    std::array<double, 3> templateReferencePoint{ 0.0, 0.0, 0.0 };
    std::vector<TargetInstanceConfig> instances;

    int globalRwgId(int targetId, int localRwgId) const {
        return targetId * rwgsPerTarget + localRwgId;
    }
};

struct IdenticalTargetAssembly {
    std::vector<Point> points;
    std::vector<FaceElement> faces;
    std::vector<RWGBase> rwgs;
    IdenticalTargetLayout layout;
};
```

全局 RWG 编号固定为：

```text
globalRwgId(a, i) = a * N0 + i
```

其中 `N0` 为模板 RWG 数。这样阻抗块、RHS、电流和 RCS 均可在 `(targetId, localRwgId)` 与全局编号间 O(1) 映射。

## 6. 用一份 NAS 构造复合几何

### 6.1 构造流程

```text
1. readNasData() 读取 templatePoints/templateFaces
2. GenerateRWG() 生成 templateRwgs
3. 从模板指针关系提取纯索引拓扑
4. 按目标数一次性 resize 复合 points/faces/rwgs
5. 复制点并施加 initialPosition
6. 复制面元并绑定到复合 points
7. 复制 RWG 并绑定到复合 faces/points
8. 设置连续 global rwgid
9. 完整性检查后构造 OCTree
```

### 6.2 指针安全

`FaceElement` 和 `RWGBase` 含有原始指针：

```text
FaceElement -> Point*
RWGBase     -> FaceElement* / Point*
OCTree Node -> RWGBase*
```

不能浅拷贝后只修改坐标，也不能在绑定指针后继续让相关 `vector` 扩容。实现要求：

- 先计算最终大小并一次性 `resize`，之后禁止导致重分配的 `push_back`。
- 用模板点索引、面索引、局部 RWG 索引重建全部指针。
- `IdenticalTargetAssembly` 的生命周期必须覆盖 OCTree 和 MLFMM。
- 禁止拷贝该对象；如允许移动，需验证缓冲区所有权和内部指针。
- 构造结束后检查所有指针均属于复合容器地址范围。

### 6.3 ID 规则

NAS 原始 ID 可能不连续，不能假定 `id == vector index + 1`。应先建立 `templatePointId -> localPointIndex`，再生成连续唯一 ID：

```text
globalPointId = targetId * pointsPerTarget + localPointIndex + 1
globalFaceId  = targetId * facesPerTarget  + localFaceIndex  + 1
globalRwgId   = targetId * rwgsPerTarget   + localRwgId
```

同时检查整数溢出。

## 7. 多目标阻抗算子

### 7.1 正确性基线

复合 `rwgs` 直接传入现有 OCTree 和 MLFMM：

```cpp
OCTree octree(assembly.rwgs, wave.wavelength(), inputFile);

MLFMM mlfmm(
    cfg, selectIntegralEqu, selectMatrixSolver,
    selectMono_Dual, polarization,
    octree.octreeNodes_, octree.octreeNodesDRvec_,
    assembly.rwgs, octree.maxLevel_,
    gausspoint, E0, wave,
    &assembly.layout);
```

全局 OCTree 会自动处理同一实例内部、不同实例之间的近区和远区作用，保留完整互耦。第一阶段不允许为追求复用而跳过任何 `a != b` 块。

### 7.2 目标内近场积分复用

纯平移同构目标满足：

```text
Z_(a,i),(a,j) = Z_template(i,j)
```

建议在 `Zmartix.h` 对同实例 near 元素采用确定性的两遍算法：

```text
第一遍：遍历全局 near pair，收集所有同实例唯一键 (localField, localSource)
第二步：并行计算每个唯一模板键的 EFIE/MFIE/CFIE 积分一次
第二遍：同实例条目从缓存复制；跨实例条目按现有函数直接计算
```

不要只复制目标 0 的 `Z_near` 行结构。不同实例相对全局 OCTree 网格的落点可能不同，所以哪些局部 RWG 对被划为近区也可能不同。缓存键集合必须取所有实例近区局部对的并集。

### 7.3 跨目标近场复用

不同实例满足：

```text
Z_(a,i),(b,j) = F(i, j, d_a - d_b)
```

只有多组目标对具有完全相同相对位移时才可复用，缓存键可为 `(relativePositionClass, localField, localSource)`。该优化留到第二版；首版直接计算跨实例近场积分，避免浮点量化和位移方向错误。

### 7.4 远区预计算

现有 `Vsmi/Vfmj` 依赖 `r - leafCenter`。目标位移若不是叶盒长度的整数倍，不同实例相对叶中心的局部坐标通常不同，不能直接复制叶层核。

第一版应：

- 正常计算复合几何的全局 `Vsmi/Vfmj`；
- 继续使用现有按节点位移向量复用的 `TF`；
- 不复制实例间叶层 `Vsmi/Vfmj`。

若性能分析证明该部分成为瓶颈，再考虑“同构目标子树 + 编队上层树”的两级森林结构，不与首版同时实现。

## 8. 多目标 RHS 与 ACA

### 8.1 平移相位关系

当前 `Vm.h` 使用：

```text
exp(-j * k * k_hat · r)
```

因此纯平移实例严格满足：

```text
b_(a,i)(t) = b_template_i(t)
              * exp(-j * k * k_hat(t) · delta_a)
```

这里 `b_template` 直接在 NAS 原坐标上计算，`delta_a = initialPosition[a] - templateReferencePoint`。该关系适用于当前 EFIE 和 CFIE RHS，因为平移不会改变 RWG 方向、面法向或极化投影。

### 8.2 `computeColumn` 优化

ACA 请求一列时只对模板 `N0` 条 RWG 做高斯积分，再用实例相位扩展为 `M*N0`：

```cpp
for (int i = 0; i < N0; ++i) {
    templateRhs[i] = computeTemplateVEntry(i, timeIndex);
}

for (int a = 0; a < M; ++a) {
    const Complex phase = exp(
        -J * k * dot(kHat[timeIndex], instanceTranslation[a]));
    for (int i = 0; i < N0; ++i) {
        rhs[a * N0 + i] = phase * templateRhs[i];
    }
}
```

高斯积分工作量由约 `O(M*N0)` 降为 `O(N0)`，额外只需复数相位乘法。

### 8.3 `computeEntry` 优化

全局行号 `g` 映射为：

```text
targetId = g / N0
localRwg = g % N0
```

返回：

```text
computeTemplateVEntry(localRwg, t) * phase(targetId, t)
```

`RCSUtils::computeMonoStatic_ACARHS_Moving` 需要获得 `IdenticalTargetLayout` 和模板 RWG，不能继续对复合重复几何逐实例做相同积分。

### 8.4 ACA 秩风险

多目标 RHS 不一定保持单目标的低秩。间距越大、运动角域越宽，实例平移相位变化越快，ACA 秩可能升高。日志必须输出：

```text
target_count
rwgs_per_target
total_rwgs
formation_aperture_m
formation_aperture_lambda
aca_rank
aca_rank / min(total_rwgs, time_count)
RHS 抽样重构误差
```

不能假定目标数增加后 rank 仍接近单目标。

## 9. 电流与 RCS

总未知量变为：

```text
N_total = targetCount * rwgsPerTarget
```

ACA 模态求解和逐时刻电流合成公式不变。联合求得的 `I(t)` 已包含互耦，不能用独立单目标电流简单拼接。

现有 `calculateRCS` 遍历复合 RWG 坐标，能够直接得到正确相干总场。进一步优化时可新增 `calculateRCSIdenticalTargets`，对观察方向 `r_hat` 复用模板远场核：

```text
F_total(t)
= sum_a exp(+j * k * r_hat(t) · delta_a)
         * sum_i F_template_i(t) * I_(a,i)(t)
```

这与复合坐标逐 RWG 积分等价，但高斯积分次数由 `M*N0` 降为 `N0`。主输出必须是全编队的相干总 RCS；各目标 dB RCS 不能直接相加。

## 10. Moving 路径语义

复合几何固定在编队本体系，实例位置就是 `initialPosition[a]`。现有 `MovingTarget::computeMovingAngleSamples` 继续把雷达方向变换到编队本体系，得到 `thetaBodyDeg(t), phiBodyDeg(t)`，随后对整个复合几何计算 RHS 和 RCS。

多目标模式下，路径行：

```text
time_s x y z qw qx qy qz
```

表示编队参考点和整体姿态，而不是某个实例的单独轨迹。建议额外导出：

```text
*_MULTI_TARGET_LAYOUT_*.txt
```

表头：

```text
target_id  initial_x  initial_y  initial_z  rwg_begin  rwg_end
```

## 11. 配置与接口

### 11.1 `RCSExportConfig.h`

建议增加：

```cpp
bool useIdenticalMultiTarget{ false };
int identicalTargetCount{ 1 };
bool reuseSelfNearIntegrals{ true };
bool reuseTemplateRHS{ true };
bool reuseTemplateFarField{ true };
```

实例位置保存在 `IdenticalTargetLayout`，不塞入 RCS 输出配置。校验规则：

```text
useIdenticalMultiTarget == true
→ selectAlgo == MLFMA
→ material == PEC
→ integralEqu == EFIE 或 CFIE
→ monoDual == Mono
→ targetCount >= 2
→ useACARHSMoving 或 useAccurateMLFMAMoving
```

首版同时接入 `Accurate-MLFMA-Moving`，用作 ACA 的逐时刻基准。

### 11.2 `MLFMM` 构造函数

增加可空参数：

```cpp
const IdenticalTargetLayout* identicalLayout = nullptr
```

求解器只读使用该布局完成 near 缓存、全局/局部 RWG 映射、RHS 相位扩展、RCS 模板复用和日志输出。单目标传 `nullptr`，旧路径不变。

### 11.3 `main.cpp`

改造入口：

```text
读取一份 NAS并生成 templateRwgs
    ↓
if useIdenticalMultiTarget:
    buildIdenticalTargetAssembly(...)
    activeRwgs = assembly.rwgs
    activeLayout = &assembly.layout
else:
    activeRwgs = templateRwgs
    activeLayout = nullptr
    ↓
OCTree(activeRwgs) → MLFMM(..., activeLayout)
```

不能通过给 `nasFiles` 重复添加同一路径来模拟此功能；那会生成多个独立任务，而不是一个联合阻抗系统。

### 11.4 批处理

用户要求在程序中给出位置，所以第一版先支持 `main.cpp` 单次模式。验证通过后再给 `BatchRunner.h` 增加可选解析：

```text
identical_target_positions = 0,0,0; 5,0,0; -5,0,0
```

### 11.5 Visual Studio 工程

新增源文件后同步更新：

```text
Code/MLFMA.vcxproj
Code/MLFMA.vcxproj.filters
```

## 12. 几何与 OCTree 检查

### 12.1 重叠检查

至少用实例 AABB 检查重叠。重叠时直接报错。非常接近但未重叠时，输出最小 AABB 间距及其波长比，提示 near 存储和迭代次数可能显著增加。

### 12.2 根盒层数问题

当前 `OCTree::computeBoundingBox` 把 `maxLevel_` 限制在 15，随后用 `rootLen = leafTarget * 2^maxLevel` 重新计算根盒。若多目标跨度要求超过 15 层，重算的根盒可能小于真实包围盒，叶索引又被 clamp，造成错误分箱。

多目标实现前必须二选一：

1. 所需层数超过上限时明确报错并输出最大场景尺寸；或
2. 在 Morton 编码安全范围内提高层数上限并验证位扩展。

建树前还要验证每个 RWG 中心都位于根盒内，禁止用索引 clamp 掩盖根盒不足。

### 12.3 完整性检查

```text
rwgs.size() == targetCount * rwgsPerTarget
rwgs[g].rwgid == g
所有点/面/RWG 指针归属正确
各实例面积、边长、法向与模板一致
各实例对应坐标差严格等于 initialPosition 差
```

## 13. 分阶段实施

### 阶段 A：正确的联合求解基线

1. 新增复合几何结构和构造器。
2. 在 `main.cpp` 增加目标位置数组。
3. 一份 NAS 生成全局连续 RWG。
4. 建立一个全局 OCTree。
5. 暂不做阻抗复用，沿用现有 `Zmartix` 和 MLFMA。
6. 让 `Accurate-MLFMA-Moving` 处理复合几何。

完成标准：结果与人工合并多目标 NAS 等价。

### 阶段 B：接入 ACA-RHS-Moving

1. 用 `row=M*N0` 接入现有 ACA。
2. 先使用全局 `computeV/computeVEntry` 验证。
3. 输出多目标 rank 和 RHS 重构误差。
4. 与 Accurate 路径逐时刻比较。

### 阶段 C：同构复用

按风险从低到高加入：

1. 模板 RHS 积分 + 实例入射相位。
2. 模板远场积分 + 实例观察相位。
3. 目标内 `Z_near` 精确积分缓存。
4. 规则阵列的跨目标相对位移缓存（可选）。

每项优化都保留关闭开关，并与未优化复数结果比较。

### 阶段 D：性能与健壮性

1. 构造前估算内存。
2. 分阶段输出几何、OCTree、`Vsmi/Vfmj`、`Z_near`、ACA、求解和 RCS 耗时。
3. 检查不同目标数和间距下的 GMRES 收敛性。
4. 根据性能数据决定是否开发两级同构子树。

## 14. 验证方案

### 14.1 单目标退化

`targetCount=1`、位置为零时，应与当前单目标的 RHS、电流和 RCS 一致。

### 14.2 双目标基准

用一份 NAS 生成两个实例，并准备人工合并的双目标 NAS。比较点坐标、面元、RWG、OCTree、RHS、电流和 RCS。

### 14.3 阻抗块

小网格用 MoM 或单位向量探测 MLFMA 算子，验证：

```text
Z00 与 Z11 的同局部索引元素一致
Z01/Z10 非零并符合互易性预期
开启/关闭 self-near 缓存结果一致
```

建议目标：

```text
模板复用前后 RHS 复数相对误差 < 1e-12
模板复用前后 near Z 相对误差 < 1e-12
联合求解电流相对误差由迭代容差控制，建议 < 1e-4
非深零点 RCS 误差 < 0.1 dB
```

### 14.4 相位与阵列干涉

设置两个远间距目标，在弱互耦情况下与“单目标复场 × 阵列因子”近似比较，用于检查 RHS 负号相位和观察场正号相位。

### 14.5 刚性运动不变性

编队整体平移且姿态不变时，在远场平面波模型下，RCS 幅值和内部干涉结构应保持不变；若输出绝对复场，则允许出现整体平移相位。

### 14.6 ACA 对比

同一路径同时运行 Accurate 和 ACA-RHS，比较 RHS、电流、线性复场、RCS dB、rank、求解次数和耗时。深零点对 dB 极敏感，必须同时检查线性复场误差。

## 15. 性能预期

设模板 RWG 数 `N0`、目标数 `M`、总未知量 `N=M*N0`：

| 项目 | 普通复合计算 | 同构优化 |
|---|---:|---|
| 几何输入 | M 份 NAS | 1 份 NAS + M 个轻量实例 |
| RHS 高斯积分/列 | `O(M*N0)` | `O(N0)` + 相位乘法 |
| RCS 高斯积分/时刻 | `O(M*N0)` | `O(N0)` + 电流/相位合成 |
| 目标内 near 积分 | 约随 M 增长 | 唯一局部 near pair 只算一次 |
| 跨目标 near 积分 | 取决于间距 | 首版仍直接计算 |
| MLFMA matvec | 约 `O(N log N)` | 保持联合算子和完整互耦 |
| ACA 求解 | rank 次 | 不变，但 rank 可能随孔径增加 |

几何和未知量仍按 `M` 线性增长，因为每个实例必须具有独立电流，才能表示互耦造成的电流差异。

## 16. 风险与措施

| 风险 | 影响 | 措施 |
|---|---|---|
| 独立运动误用固定 Z | 物理模型错误 | 首版只允许刚性编队 |
| 忽略互耦 | 多目标结果错误 | 一个全局 OCTree 和联合求解器 |
| vector 扩容导致指针失效 | 崩溃或隐蔽错误 | 先 resize 再绑定，增加地址检查 |
| NAS ID 不连续 | 错误连接 | 显式 ID 到局部索引映射 |
| OCTree 根盒不足 | 错误分箱 | 修复层数上限并验证范围 |
| 错误复制 `Vsmi/Vfmj` | 叶中心相位错误 | 首版正常计算全局叶层核 |
| 平移相位符号错误 | 干涉条纹错误 | 与直接积分、阵列因子双重对比 |
| 编队孔径使 ACA rank 上升 | 加速比下降 | 输出孔径、rank、误差，不强截断 |
| 目标过近或重叠 | near 和迭代代价激增 | AABB/间距检查，重叠报错 |

## 17. 预计修改文件

| 文件 | 内容 |
|---|---|
| `Code/IdenticalTargetAssembly.h/.cpp` | 实例、布局、复合几何、指针重建和检查 |
| `Code/main.cpp` | 目标位置数组和多目标入口 |
| `Code/RCSExportConfig.h` | 多目标开关、数量、输出标签和校验 |
| `Code/MLFMM.h/.cpp` | 接收布局并传入 near/RHS/RCS |
| `Code/Zmartix.h` | 同实例 near 积分两遍缓存 |
| `Code/RCS.h` | 模板 RHS、模板远场和布局输出 |
| `Code/OCTree.cpp` | 大场景根盒层数与越界检查 |
| `Code/BatchRunner.h` | 后续增加位置列表批处理 |
| `Code/MLFMA.vcxproj(.filters)` | 注册新增文件 |

`NASReader` 和 `RWGGenerator` 的单目标行为不做破坏性修改，多目标复制由独立模块承担。

## 18. 审核确认项

1. 多目标运动按刚性编队处理，相对位置始终不变。
2. 所有目标使用同一 NAS、相同尺度、材料和朝向，首版只给平移位置。
3. `initialPosition` 表示 NAS 模板参考点在编队坐标系中的位置，默认参考点为 NAS 原点。
4. 使用一个全局 MLFMA 联合系统并保留全部互耦。
5. 首版只支持 `PEC + Mono + EFIE/CFIE + MLFMA`。
6. 运动路径描述编队参考点，而不是单个实例。
7. 先实现与人工合并 NAS 等价的基线，再开启 RHS、RCS、near 复用。
8. 首版不复制不同实例的 `Vsmi/Vfmj`。
9. 主输出是全体目标的相干总 RCS，不直接相加各目标 dB。
10. 独立轨迹、不同实例姿态、球面波雷达和时变阻抗算子留待后续版本。

审核通过后，建议严格按阶段 A、B、C、D 实施，使几何复制、互耦、ACA 和相位复用的误差可以分别验证。
