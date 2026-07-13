# MLFMM 插值函数与 Fortran 版本对比分析

## 1. 对比对象

本次对比的 C++ 代码位于：

```txt
Code/AlgoMLFMM.h
```

重点函数为：

```cpp
AlgoMLFMM::find_interp_pos(...)
AlgoMLFMM::computeInterpolationMatrix(...)
lagrangeWeight(...)
```

Fortran 参考代码位于：

```txt
D:\MyCode\Integral_Solver_V3.2_OMP\Integral_Solver_V3.2_OMP\Integral Solver Kernel\MLFMA\MLFMA.f90
```

重点过程为：

```fortran
set_interp_CFIE
set_interp_PMCHW
find_interp_pos
weight_cal_function
interpolation
anterpolation
```

## 2. C++ 当前实现

C++ 中的插值点查找函数为：

```cpp
inline std::array<int, 3> find_interp_pos(
    double target,
    const std::vector<double>& arr,
    bool periodic = true)
```

其主要逻辑是：

1. 根据 `arr[1] - arr[0]` 判断数组升序或降序。
2. 用二分搜索找到离 `target` 最近的点 `closest`。
3. 返回三个插值点：

   ```txt
   prev, curr, next
   ```

4. 如果 `periodic == true`，边界点会回绕：

   ```txt
   closest = 0      -> prev = n - 1
   closest = n - 1  -> next = 0
   ```

5. 如果 `periodic == false`，当前实现会在边界返回重复点：

   ```cpp
   closest == 0      -> {0, 0, 1}
   closest == n - 1  -> {n - 2, n - 1, n - 1}
   ```

这会导致三点 Lagrange 插值分母出现 0，因此当前 `periodic == false` 不能直接安全用于二次插值。

在 `computeInterpolationMatrix()` 中，当前代码为：

```cpp
auto theta_idx = find_interp_pos(theta0, theta_level[revS], true);
auto phi_idx   = find_interp_pos(phi0,   phi_level[revS],   true);
```

也就是说，当前 C++ 代码把 `theta` 和 `phi` 都当成周期变量处理。

## 3. Fortran 参考实现

Fortran 的 `find_interp_pos` 接口为：

```fortran
function find_interp_pos(value_true,mydata) result(pos)
```

它没有显式的 `periodic` 参数。主要逻辑为：

1. 根据：

   ```fortran
   delta = mydata(2) - mydata(1)
   ```

   判断升序或降序。

2. 用二分搜索找到离 `value_true` 最近的 `midloc`。

3. 返回三个插值位置：

   ```fortran
   pos(1), pos(2), pos(3)
   ```

4. 在边界处直接回绕：

   升序时：

   ```fortran
   if(midloc == 1) then
       pos(1) = num
       pos(3) = midloc + 1
   else if(midloc == num) then
       pos(1) = midloc - 1
       pos(3) = 1
   end if
   ```

   降序时：

   ```fortran
   if(midloc == 1) then
       pos(1) = midloc + 1
       pos(3) = num
   else if(midloc == num) then
       pos(1) = 1
       pos(3) = midloc - 1
   end if
   ```

因此，Fortran 版本的 `find_interp_pos` 默认就是周期边界版本。

Fortran 权重函数为：

```fortran
function weight_cal_function(x,x_guess) result(weight)
```

其公式是标准三点 Lagrange 权重：

```fortran
tmp = tmp * (x - x_guess(j)) / (x_guess(i) - x_guess(j))
```

C++ 中的 `lagrangeWeight(...)` 与该公式等价。

## 4. set_interp_PMCHW 对 PMCHWT 的处理方式

Fortran 的 `set_interp_PMCHW` 分别对外部空间和介质内部空间构造插值表。

外部空间：

```fortran
interp_pos_phai(:)  = find_interp_pos(phai_real,  this%octree(i+1)%sph_ang_phai)
interp_pos_theta(:) = find_interp_pos(theta_real, this%octree(i+1)%sph_ang_theta)
```

介质内部空间：

```fortran
interp_pos_phai(:)  = find_interp_pos(phai_real,  this%octree(i+1)%sph_ang_phai2)
interp_pos_theta(:) = find_interp_pos(theta_real, this%octree(i+1)%sph_ang_theta2)
```

也就是说，Fortran 版本对 `theta` 和 `phai` 都调用同一个周期式 `find_interp_pos`。

C++ 当前实现也做了等价处理：

```cpp
auto theta_idx = find_interp_pos(theta0, theta_level[revS], true);
auto phi_idx   = find_interp_pos(phi0,   phi_level[revS],   true);
```

从“是否忠实移植 Fortran”角度看，C++ 当前写法与 Fortran 是一致的。

## 5. 关键差异

### 5.1 索引基准不同

Fortran 是 1-based：

```fortran
pos(1), pos(2), pos(3)
```

C++ 是 0-based：

```cpp
return { prev, curr, next };
```

C++ 已经在 `InterpData` 中按 0-based 方式计算：

```cpp
tp[k++] = tIdx[i] * phiCount + pIdx[j];
```

这一点没有明显错误。

### 5.2 C++ 增加了 periodic 参数

Fortran 没有 `periodic` 参数，默认所有方向都周期回绕。

C++ 增加了：

```cpp
bool periodic = true
```

但当前 `computeInterpolationMatrix()` 中仍然给 `theta` 和 `phi` 都传 `true`。因此实际行为仍接近 Fortran。

### 5.3 C++ 的非周期分支不能直接用于 Lagrange

C++ 的 `periodic == false` 分支在边界会产生重复点：

```cpp
closest == 0      -> {0, 0, 1}
closest == n - 1  -> {n - 2, n - 1, n - 1}
```

如果后续调用：

```cpp
lagrangeWeight(theta0, theta_vals[0], theta_vals[1], theta_vals[2])
```

就可能出现：

```txt
theta_vals[0] == theta_vals[1]
```

从而导致除零。因此如果要把 `theta` 改成非周期插值，必须同时修复非周期边界策略。

### 5.4 phi 周期角度没有展开

Fortran 和 C++ 都有一个共同风险：`phai/phi` 在索引上做了周期回绕，但权重计算时没有把角度值展开到同一个连续区间。

例如：

```txt
phi0 ≈ 0
phi_vals = {接近 2π, 接近 0, 接近 Δφ}
```

直接做 Lagrange：

```txt
L(phi0; 2π - Δφ, Δφ/2, 3Δφ/2)
```

会把本来相邻的周期点当成距离很远的点。

这一点 Fortran 原版也存在，但如果旧代码在测试范围内表现可接受，可能是因为采样点避开了最敏感位置，或者误差尚未在低频场景中放大。

## 6. 对 MLFMM 不收敛问题的解释

如果目标是解释为什么 C++ 的 MLFMM 在某些介质球、较高频率或更细网格下不收敛，那么需要区分两个层次。

### 6.1 从移植一致性看

C++ 当前 `find_interp_pos(..., true)` 的行为与 Fortran 的 `find_interp_pos` 高度一致：

- 都支持升序和降序。
- 都找最近点作为中点。
- 都在边界使用周期回绕。
- 都使用三点 Lagrange 权重。

因此，如果只问“C++ 是否按 Fortran 写法移植”，答案是：基本是的。

### 6.2 从数学合理性看

`phi` 是周期变量，周期回绕合理。

但 `theta` 的范围是 `[0, π]`，不是普通周期变量。把 `theta = 0` 附近与 `theta = π` 附近直接作为相邻插值点，在球坐标角度上并不自然。虽然球面极点附近存在方向退化，但直接使用一维 `theta` 周期 Lagrange 插值容易引入误差。

因此，C++ 忠实继承了 Fortran 的写法，也继承了它的潜在风险。这个风险在以下情况下更容易暴露：

- 频率升高，电尺寸变大。
- 介质内部波数 `k2` 更大。
- 多层层数增加，`m2m/l2l` 插值次数增加。
- PMCHWT 系统条件数较差。
- GMRES 或 CGS 对矩阵向量乘法误差更敏感。

这与 `Sphere_Die_8e8.nas` 能得到结果，而 `Sphere_Die_1e9.nas` 迭代困难的现象是吻合的。

## 7. 建议修改方案

### 7.1 保守方案：保持 Fortran 一致

如果你的目标是严格复现 Fortran 代码行为，C++ 当前 `find_interp_pos` 不需要大改。

但建议补充注释，明确：

```cpp
// This follows the original Fortran implementation:
// both theta and phi use periodic neighbor selection.
```

这样以后排查时不会误以为 `theta` 周期处理是无意错误。

### 7.2 推荐方案：区分 theta 与 phi

如果你的目标是提升 C++ MLFMM 在高频/介质 PMCHWT 下的稳定性，建议改成：

```cpp
auto theta_idx = find_interp_pos(theta0, theta_level[revS], false);
auto phi_idx   = find_interp_pos(phi0,   phi_level[revS],   true);
```

同时修复 `find_interp_pos(..., false)` 的边界逻辑，确保返回三个不重复点：

```cpp
if (!periodic) {
    if (closest == 0) {
        return { 0, 1, 2 };
    }
    if (closest == n - 1) {
        return { n - 3, n - 2, n - 1 };
    }
    return { closest - 1, closest, closest + 1 };
}
```

### 7.3 phi 权重建议做角度展开

对 `phi`，索引可以周期，但用于权重计算的角度值建议展开：

```cpp
double phi_vals[3] = {
    phi_level[revS][phi_idx[0]],
    phi_level[revS][phi_idx[1]],
    phi_level[revS][phi_idx[2]]
};

for (double& v : phi_vals) {
    if (v - phi0 > Pi) {
        v -= 2.0 * Pi;
    }
    if (phi0 - v > Pi) {
        v += 2.0 * Pi;
    }
}
```

这样跨 `0/2π` 边界时，Lagrange 插值看到的是连续角度。

## 8. 建议验证方法

不要只通过 RCS 判断插值是否正确，建议直接比较矩阵向量乘法。

对同一个小模型、同一个随机向量 `x`，分别计算：

```txt
y_fmm   = A_fmm   * x
y_mlfmm = A_mlfmm * x
```

比较：

```txt
||y_mlfmm - y_fmm|| / ||y_fmm||
```

建议分三组测试：

1. 当前 C++ 写法：`theta periodic = true`，`phi periodic = true`。
2. 修改后：`theta periodic = false`，`phi periodic = true`。
3. 修改后并加入 `phi` 角度展开。

如果第 2 或第 3 组显著降低 `MLFMM` 与 `FMM` 的 `Ax` 差异，那么可以确认插值边界处理是导致高频介质 PMCHWT 不收敛的重要原因。

## 9. 结论

从代码对比看，C++ 的 `find_interp_pos` 当前实现基本是 Fortran `find_interp_pos` 的 0-based 版本，并且当前 `computeInterpolationMatrix()` 也沿用了 Fortran 对 `theta/phai` 都周期回绕的策略。

但是，从 MLFMM 数值稳定性看，`theta` 周期插值和 `phi` 跨边界角度未展开都有潜在风险。它们在低频或较粗场景下可能不明显，但在 `Sphere_Die_1e9.nas` 这类更高频、PMCHWT 介质问题中更容易放大为矩阵向量乘法误差，最终表现为 MLFMM 迭代不收敛或残差停滞。

建议下一步优先做 `FMM Ax` 与 `MLFMM Ax` 的直接对比，并测试 `theta` 非周期插值加 `phi` 角度展开后的误差变化。
