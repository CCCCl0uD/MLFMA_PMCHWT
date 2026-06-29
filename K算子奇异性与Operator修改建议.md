# K 算子奇异性判断与 `Operator.h` 修改建议

## 1. 结论

对当前代码所实现的 Galerkin 矩阵元

\[
K_{mn}(k)=\int_{S_m}\boldsymbol\Lambda_m(\mathbf r)\cdot
\int_{S_n}\nabla_{\mathbf r}G_k(\mathbf r,\mathbf r')
\times\boldsymbol\Lambda_n(\mathbf r')\,\mathrm dS'\,\mathrm dS,
\]

不能笼统地说“`K` 算子没有奇异问题”。更准确的结论是：

1. **核函数本身有强奇异性。** 当 \(R=|\mathbf r-\mathbf r'|\to0\) 时，
   \(\nabla G_k=O(R^{-2})\)。因此不能把任意相交或相邻三角形都当作普通光滑积分。
2. **同一平面三角形的本公式矩阵元逐点为零。** 对 `rwgFieldTriId == rwgSourceTriId`，
   \(\boldsymbol\Lambda_m\)、\(\boldsymbol\Lambda_n\) 和
   \(\mathbf R=\mathbf r-\mathbf r'\) 都位于同一切平面，故
   \[
   \boldsymbol\Lambda_m\cdot(\mathbf R\times\boldsymbol\Lambda_n)=0.
   \]
   这个结论与 \(R\) 的大小无关，所以同三角形项应直接置零，不需要当前的 `singK` 解析补偿。
3. **共边、共点和距离很近但不相同的三角形仍有奇异或近奇异数值问题。** 它们不满足上述“全部向量共面”的逐点抵消。积分可能以反常积分意义存在，但普通固定阶高斯积分不一定准确。
4. **PMCHWT 实际使用的是 \(K(k_1)-K(k_2)\)。** 两项的波数无关主奇异部分完全相消，差分核只剩有界的向量核。直接计算差分比先分别计算两个大数再相减更稳定，也更省时。

因此，若“没有奇异问题”专指当前公式的**同一个三角形自项**，说法正确；若指单独的 \(K(k)\) 在所有几何关系下都无奇异性，则说法不正确。

> 本结论只适用于文首这个未经法向量旋转测试的公式。如果以后实现的是
> \(\hat{\mathbf n}\times K\)、MFIE 极限算子，或显式包含 \(\pm\tfrac12 I\) 跳跃项的定义，不能直接把自项置零。

## 2. 与当前代码的对应关系

`Code/Operator.h` 中：

- 第 975 行的标量系数
  ```cpp
  exp(-J * k_ * R) * (1.0 + J * k_ * R) / (R * R * R)
  ```
  与 `Rvec` 相乘后就是（差一个由整体 `1/(4*pi)` 处理的常数和符号约定）\(\nabla G_k\)。
- 第 980–982 行计算
  ```cpp
  dot(pomocvF, cross(Rvec, pomocvS))
  ```
  正是标量三重积
  \(\boldsymbol\Lambda_m\cdot(\mathbf R\times\boldsymbol\Lambda_n)\)。
- 第 987–1036 行对同一三角形先拆除 `1/R^3` 与 `k^2/(2R)`，再调用 `singK` 加回。
- 第 324–706 行的 `singK` 就是上述解析补偿实现。

问题在于：同一三角形上标量三重积理论上恒为零。当前代码仍执行双重高斯积分和一大段解析补偿，不仅没有必要，还可能因浮点投影、法向量误差和多个大数相消而产生本不应存在的小非零自项。

## 3. 数学推导

采用代码中的格林函数约定（整体 \(1/(4\pi)\) 在函数末尾乘入）：

\[
G_k(R)=\frac{e^{-\mathrm j kR}}{R},\qquad
\nabla_{\mathbf r}G_k
=-\frac{e^{-\mathrm j kR}(1+\mathrm j kR)}{R^3}\mathbf R.
\]

符号取决于使用 \(\nabla_{\mathbf r}\) 还是 \(\nabla_{\mathbf r'}\)，但不影响奇异阶数及以下判断。小量展开为

\[
\frac{e^{-\mathrm j kR}(1+\mathrm j kR)}{R^3}
=\frac1{R^3}+\frac{k^2}{2R}-\mathrm j\frac{k^3}{3}+O(R).
\]

乘上 \(\mathbf R\) 后，单个 \(K(k)\) 的主核为 \(O(R^{-2})\)。这说明核并非光滑。

但在同一平面三角形上，RWG 基函数是切向量，并且
\(\mathbf R\) 也是切向量。于是 \(\mathbf R\times\boldsymbol\Lambda_n\)
垂直于三角形，而 \(\boldsymbol\Lambda_m\) 位于三角形内，两者点积严格为零：

\[
\boldsymbol\Lambda_m\cdot
(\mathbf R\times\boldsymbol\Lambda_n)=0.
\]

对于 PMCHWT 差分，展开式给出

\[
\frac{e^{-\mathrm j k_1R}(1+\mathrm j k_1R)
-e^{-\mathrm j k_2R}(1+\mathrm j k_2R)}{R^3}
=\frac{k_1^2-k_2^2}{2R}
-\mathrm j\frac{k_1^3-k_2^3}{3}+O(R).
\]

再乘 \(\mathbf R\) 后为 \(O(1)\)，所以 `K1 - K2` 的向量核有界；不过方向极限一般依赖接近路径，仍应避免在 `R == 0` 处直接除法。

## 4. 修改方案 A：最小且安全的修改

这是最适合先实施和验证的方案。保留 `K_operator(k)` 接口及所有调用，仅修正同三角形分支。

### 4.1 用“自三角形贡献为零”替换当前同面分支

在 `K_operator` 中，将当前第 987–1037 行整个 `else` 分支替换为：

```cpp
else {
    // 对当前 Galerkin K 公式，同一平面三角形上
    // Lambda_m、R = r-r'、Lambda_n 均为切向量，因此
    // Lambda_m · (R x Lambda_n) == 0。
    // 自三角形贡献理论上严格为零；不做奇异拆分与补偿。
    alok1 = std::complex<double>(0.0, 0.0);
}
```

实际上 `alok1` 在进入分支前已经初始化为零，因此也可以写成：

```cpp
else {
    // Same-triangle Galerkin K contribution is identically zero.
}
```

建议保留显式赋零和数学原因注释，以防以后有人改变初始化位置。

### 4.2 删除 `singK`

确认工程内没有其他调用后，可删除 `Operator.h` 第 324–706 行的完整 `singK` 函数。当前仓库中它只在 `K_operator` 的同三角形分支被调用。

如果希望分两步降低改动风险，可以先只停止调用并保留函数；数值回归通过后再删除死代码。

### 4.3 给非同面分支增加退化网格保护

不同 `faceId` 正常情况下不应在高斯内点上出现 `R == 0`。若网格存在重复重合面，当前第 975 行会直接除零。建议在计算核之前加入诊断：

```cpp
if (R <= eps0) {
    // 不同 faceId 却出现重合积分点，通常意味着重复面或退化网格。
    // 不能静默把一个真正的奇异项设为零。
    throw std::runtime_error(
        "K_operator: coincident quadrature points on different faces");
}
```

同时需要在文件头加入：

```cpp
#include <stdexcept>
```

若生产版本不允许抛异常，至少应使用断言并在网格预处理阶段检查重复三角形。不要简单采用 `if (R < eps0) G = 0`，那会掩盖几何错误。

### 4.4 共边/共点与近邻面

`rwgFieldTriId != rwgSourceTriId` 并不等价于“光滑远场”。对于共边、共点和极近三角形，建议后续按几何关系选择：

- 共边三角形：Duffy 变换、奇异性消去或专用邻接三角形积分；
- 共点三角形：顶点型 Duffy 变换或自适应细分；
- 不相交但很近：自适应细分或提高积分阶数；
- 普通分离三角形：保留当前双重高斯积分。

这部分属于精度增强，不是把同三角形分支删掉后才产生的新问题；当前实现本来就存在这项风险。

## 5. 修改方案 B：PMCHWT 推荐方案——直接计算 `K(k1)-K(k2)`

`Code/Zmartix.h` 第 514–519 行先计算 `ZK1`、`ZK2`，随后只使用二者之差。推荐新增一个专门的差分函数，避免重复遍历和主奇异项的数值相消。

### 5.1 新接口

```cpp
inline void K_difference_operator(
    const RWGBase* rwgField,
    const RWGBase* rwgSource,
    const std::complex<double> k1,
    const std::complex<double> k2,
    const gaussPoints& gp,
    std::complex<double>& Zij_near_Kdiff);
```

其几何、RWG 正负三角形和最终系数与当前 `K_operator` 完全相同，只把核系数改为

```cpp
const std::complex<double> x1 = J * k1 * R;
const std::complex<double> x2 = J * k2 * R;

const std::complex<double> numerator =
    std::exp(-J * k1 * R) * (1.0 + x1) -
    std::exp(-J * k2 * R) * (1.0 + x2);

const std::complex<double> kernelDiff = numerator / (R * R * R);
```

然后继续使用：

```cpp
cross(&CROSS1[0], &Rvec[0], &pomocvS[0]);
alok1 += wF * wS * kernelDiff * dot(&pomocvF[0], &CROSS1[0]);
```

同一三角形仍直接置零。

### 5.2 小距离下避免灾难性消减

上述 `numerator` 在很小的 \(R\) 下是两个接近 1 的复数之差，可能损失有效数字。建议根据无量纲量
\(\max(|k_1R|,|k_2R|)\) 切换到级数：

```cpp
inline std::complex<double> KKernelDifference(
    const std::complex<double> k1,
    const std::complex<double> k2,
    const double R)
{
    const double scale = std::max(std::abs(k1 * R), std::abs(k2 * R));

    if (R <= eps0) {
        // 完整被积式还会乘 Rvec；在孤立的重合采样点取零。
        return std::complex<double>(0.0, 0.0);
    }

    if (scale < 1.0e-3) {
        const std::complex<double> dk2 = k1 * k1 - k2 * k2;
        const std::complex<double> dk3 = k1 * k1 * k1 - k2 * k2 * k2;
        const std::complex<double> dk4 =
            k1 * k1 * k1 * k1 - k2 * k2 * k2 * k2;

        return 0.5 * dk2 / R
            - (J / 3.0) * dk3
            - 0.125 * dk4 * R;
    }

    return (
        std::exp(-J * k1 * R) * (1.0 + J * k1 * R) -
        std::exp(-J * k2 * R) * (1.0 + J * k2 * R)
    ) / (R * R * R);
}
```

注意：函数返回的标量系数在 \(R\to0\) 时含 `1/R`，但它随后与 `Rvec` 叉乘，所以完整向量核有界。`R == 0` 只是数值积分中的单个采样点，按零处理不会改变积分值；若它来自不同编号的重合面，则仍应优先报告网格错误。

### 5.3 修改 `Zmartix.h`

把：

```cpp
std::complex<double> ZK1(0.0, 0.0), ZK2(0.0, 0.0);

OLK::K_operator(rwgField, rwgSource, wave.k1(), gausspoint, ZK1);
OLK::K_operator(rwgField, rwgSource, wave.k2(), gausspoint, ZK2);

std::complex<double> z12 = ZK1 - ZK2;
std::complex<double> z21 = -ZK1 + ZK2;
```

改为：

```cpp
std::complex<double> ZKdiff(0.0, 0.0);

OLK::K_difference_operator(
    rwgField, rwgSource,
    wave.k1(), wave.k2(),
    gausspoint, ZKdiff);

std::complex<double> z12 = ZKdiff;
std::complex<double> z21 = -ZKdiff;
```

`ZL1`、`ZL2` 和对角块的代码保持不变。

方案 B 只适用于调用方确实只需要差分的地方。如果别处需要单独的 \(K(k)\)，仍应保留单波数接口，并对共边/共点积分做适当处理。

## 6. 建议实施顺序

1. 先采用方案 A：同三角形直接置零，暂时保留但不调用 `singK`。
2. 对小网格保存修改前后的矩阵，检查所有同 `faceId` 子贡献是否从浮点小量变为严格零。
3. 比较最终 PMCHWT 矩阵和散射结果。若变化明显，应优先排查旧 `singK` 是否产生了伪自项，而不是恢复它。
4. 数值回归通过后删除 `singK` 死代码。
5. 再实施方案 B，比较 `ZKdiff` 与旧的 `ZK1-ZK2`；中远距离应在舍入误差内一致，小距离通常会更稳定。
6. 最后再引入共边、共点和近邻面的专用积分策略，并做网格加密收敛测试。

## 7. 验证建议

至少增加以下测试：

### 7.1 同三角形代数恒等式

随机选择同一三角形内部的两点，构造两个 RWG 线性向量，验证

```cpp
std::abs(dot(RWGf, cross(Rvec, RWGs))) < tolerance
```

其中 `tolerance` 应按三角形尺度归一化，而不是使用固定绝对误差。

### 7.2 自项严格为零

让 field/source 子三角形 `faceId` 相同，检查 `alok1 == 0`。同时改变高斯积分阶数，结果应保持严格为零。

### 7.3 差分核一致性

对中等距离的随机点比较：

```cpp
KKernelDifference(k1, k2, R)
```

与旧实现的两个核系数之差。距离较小时使用高精度参考值验证级数分支。

### 7.4 邻接与近邻收敛

分别选取共边、共点、近邻和普通分离三角形，提高积分阶数或细分层数，观察矩阵元是否收敛。这个测试用于判断何时必须引入 Duffy 变换或自适应积分。

### 7.5 物理量回归

在已有介质散射算例上比较：

- PMCHWT 矩阵条件数或迭代残差；
- 表面电流、磁流；
- RCS 或远场；
- 网格加密前后的收敛趋势。

## 8. 最终建议

当前最应修改的是：**删除 `K_operator` 同三角形分支中的数值积分和 `singK` 补偿，明确令该子三角形贡献为零。** 这与给定公式完全一致。

如果该代码路径只服务于 `Zmartix.h` 中的 PMCHWT 组装，则进一步推荐直接实现 `K_difference_operator(k1, k2)`。这会在离散前消掉 \(1/R^3\) 主项，减少一半 K 积分工作量，并避免 `ZK1`、`ZK2` 两个近似大数相减造成的精度损失。

但不要把这推广成“所有 K 积分都可以用普通高斯积分”：不同但共边、共点或很近的三角形仍需要专门的数值积分或至少收敛性验证。
