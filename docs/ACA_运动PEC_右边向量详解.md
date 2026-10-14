# ACA 加速矩量法：单个运动 PEC 目标右边向量的 ACA 分解与作用

> 项目: `pmchw_fortri` | 方法: PMCHWT + ACA | 目标: 运动 PEC
> 路径: `D:\CEM\MLH\ACA-MOM-finally\pmchw_fortri`

---

## 一、整体架构

### 1.1 矩量法线性系统

```
Z * x = b
```

- `Z` (N*N): 阻抗矩阵，元素 Z_mn 表示第 n 个 RWG 基函数对第 m 个测试函数的互作用
- `b` (N*1): 右边激励向量，元素 b_m 表示入射平面波在第 m 个边上的投影
- `x` (N*1): 待求电流系数向量

### 1.2 ACA 在本项目中的两种用途

| 用途 | 压缩对象 | 矩阵维度 | 核心文件 |
|------|----------|----------|----------|
| **Z-ACA** | 阻抗矩阵远场块 | 边缘*边缘 | `ACA.cpp` |
| **V-ACA** | **右边激励向量矩阵** | **边缘*角度(或时间)** | `ACA_V.cpp`, `ACA_v_mov.cpp` |

> 本文档聚焦 **V-ACA**——即如何用 ACA 对右边激励向量做低秩压缩。

### 1.3 V-ACA 的物理含义

对于单站 RCS 或运动目标，激励向量随角度/时间变化形成矩阵:

```
B(边缘索引, 角度/时间索引)
```

- **静态单站**: 维度 = N_edges * N_angles
- **运动目标**: 维度 = N_edges * N_timesteps

ACA 将该矩阵分解为:

```
B ~ U * V^T
```

| 矩阵 | 维度 | 物理含义 |
|------|------|----------|
| U | N_edges * rank | "空间模态": 每条边在不同激励模式下的权重 |
| V | N_angles * rank | "角度/时间模态": 每个角度/时刻对各模式的贡献系数 |

---

## 二、关键参数与全局常量

文件 `mymathparater.h`:

```cpp
#define COMPLEX    complex<double>
#define I          COMPLEX(0,1)
#define PI         double(3.14159265358979323846)
#define eta0       COMPLEX(119.9169832 * PI)    // 自由空间波阻抗 ~ 377 Ohm
#define eps_ACA    double(1.0e-3)              // ACA 收敛误差容限
#define Aleph      double(0.2)                  // PMCHWT 组合系数 alpha
```

**PMCHWT 组合公式**:

```
Z_total = alpha * Z_EFIE + (1-alpha) * Z_MFIE
b_total = alpha * b_E     + (1-alpha) * b_H
```

---

## 三、底层: 单条边的激励积分计算

文件 `excvecE.cpp`，函数 `assemble_exic_vector`。

### 3.1 函数签名

```cpp
void assemble_exic_vector(
    vector<COMPLEX>& C1,          // 输出: EFIE 分量 b_E
    vector<COMPLEX>& C2,          // 输出: MFIE 分量 b_H
    Mesh& mesh, Points& points,   // 网格 + 积分点
    int n_maxele,                 // 全局棱边编号
    COMPLEX k, COMPLEX eta,       // 波数, 波阻抗
    double thinc, double phinc,   // 入射俯仰角, 方位角(度)
    double pol_i,                 // 极化角(度)
    Edge_info Edges,              // 棱边拓扑
    double* trans_VecF,           // 平动向量 [dx, dy, dz]
    vector<double>& rod_rig_F,    // 3*3旋转矩阵展开为 9 元素
    int maxele);                  // 每目标边数(用于取模定位)
```

### 3.2 计算流程

```
1. 初始化平面波: k_i(波矢量), e_i(电场方向), h_i(磁场方向)

2. f_rel = n_maxele % maxele  -> 找到局部棱边索引

3. 获取该棱边共享的两个三角形 T+, T-:
   ele1 = Edges.get_triedg(f_rel)   // 三角形编号
   topF = Edges.get_topedg(f_rel)   // 对顶点局部编号

4. 对每个三角形 (f_edg_n = 0,1):
   a) 取出顶点坐标, 施加刚体变换:  p' = R * p + trans
   b) 三角形法向量也旋转:  n' = R * n
   c) 三角形面积: AR = det(Jacobian) / 2

   d) 高斯积分 (面积坐标 N1,N2; 权重 w):
      - 场点: r_g = N1*p1' + N2*p2' + (1-N1-N2)*p3'
      - RWG 基函数: f_m(r_g) = (r_g - r_vertex') / (2*AR)
      - 相位: phase = exp(-j * k * (k_i * r_g))

      b_E += w * det * (f_m * e_i) * phase
      b_H += w * det * (f_m * (n' x h_i)) * phase

   e) 符号: s = (f_edg_n==0) ? -1 : 1
      C1 += s * b_E
      C2 += s * b_H
```

### 3.3 完整代码(核心算法)

```cpp
namespace EFIE { namespace excEFIE {

void assemble_exic_vector(
    vector<COMPLEX>& C1, vector<COMPLEX>& C2,
    Mesh& mesh, vector<Trianinfo>& Triangles, Points& points,
    int n_maxele, COMPLEX k, COMPLEX eta,
    double thinc, double phinc, double pol_i, Edge_info Edges,
    double* trans_VecF, vector<double>& rod_rig_F, int maxele)
{
    // ---- 1. 获取平面波 ----
    EMSouce EM_ki;
    EM_ki.initPW(pol_i, thinc, phinc);
    double k_i[3] = {EM_ki.vW(0), EM_ki.vW(1), EM_ki.vW(2)};
    double e_i[3] = {EM_ki.vE(0), EM_ki.vE(1), EM_ki.vE(2)};
    double h_i[3] = {EM_ki.vH(0), EM_ki.vH(1), EM_ki.vH(2)};

    // ---- 2. 获取局部索引 ----
    int f_rel = n_maxele % maxele;
    int* ele1 = Edges.get_triedg(f_rel);  // 两个三角形编号
    int* topF = Edges.get_topedg(f_rel);  // 对顶点局部编号

    // ---- 3. 遍历两个三角形 ----
    for (int f_edg_n = 0; f_edg_n < 2; f_edg_n++) {
        int nn_tri = ele1[f_edg_n];
        int nn_top = topF[f_edg_n];
        int* n1 = mesh.getNOvertex(nn_tri);
        double* p11 = mesh.getCoord(n1[0]);
        double* p21 = mesh.getCoord(n1[1]);
        double* p31 = mesh.getCoord(n1[2]);

        // 施加刚体变换: 旋转 + 平动
        vector<double> p1(3), p2(3), p3(3);
        Rodvec_tra(&p1[0], p11, &rod_rig_F[0], trans_VecF);
        Rodvec_tra(&p2[0], p21, &rod_rig_F[0], trans_VecF);
        Rodvec_tra(&p3[0], p31, &rod_rig_F[0], trans_VecF);

        // 法向量旋转
        Trianinfo trianF = Triangles[nn_tri];
        vector<double> nvecF1 = trianF.getnorm();
        vector<double> nvecF(3);
        Rodvec(&nvecF[0], &nvecF1[0], &rod_rig_F[0]);

        double det1 = trianF.getDeter();
        double AR1 = det1 / 2.0;

        // ---- 4. 高斯积分 ----
        const auto& PointsS = points.getPointsS();
        const auto& WeightsS = points.getWeightsS();
        COMPLEX alok1exvec(0), alok2exvec(0);

        for (unsigned int nf = 0; nf < WeightsS.size(); ++nf) {
            double N1 = PointsS[nf][0], N2 = PointsS[nf][1];
            double N0 = 1.0 - N1 - N2;
            double wf = WeightsS[nf];

            // 场点坐标
            double fielpoin[3];
            fielpoin[0] = p1[0]*N1 + p2[0]*N2 + p3[0]*N0;
            fielpoin[1] = p1[1]*N1 + p2[1]*N2 + p3[1]*N0;
            fielpoin[2] = p1[2]*N1 + p2[2]*N2 + p3[2]*N0;

            // RWG 基函数: f = (r - r_vertex) / (2*area)
            double* p111 = mesh.getCoord(n1[nn_top]);
            vector<double> p(3);
            Rodvec_tra(&p[0], p111, &rod_rig_F[0], trans_VecF);
            double pomocvF[3], RWGf[3];
            subtract(&pomocvF[0], &fielpoin[0], &p[0]);
            multconst(&RWGf[0], &pomocvF[0], 1.0/(2.0*AR1));

            // n x h_i
            double cross_ev[3];
            cross(&cross_ev[0], &nvecF[0], &h_i[0]);

            // 积分累加
            COMPLEX phase = exp(-I * k * dot(&fielpoin[0], k_i));
            alok1exvec += wf * det1 * dot(&RWGf[0], e_i) * phase;
            alok2exvec += wf * det1 * dot(&RWGf[0], &cross_ev[0]) * phase;
        }

        // 符号修正 (RWG 基函数定义)
        const int s1 = (f_edg_n == 0) ? -1 : 1;
        C1[0] += double(s1) * alok1exvec;
        C2[0] += double(s1) * alok2exvec;
    }
}

}} // namespace
```

### 3.4 两种坐标变换的区别

```cpp
// === assemble_exic_vector (静态/目标自身姿态) ===
// p' = R * p + trans (先旋转后平动)
Rodvec_tra(&p[0], coords, &rod_rig_F, trans_VecF)
// 实现: u[i] = Sum_j R[i*3+j]*v[j] + t[i]

// === assemble_exic_vector_mov (运动轨迹版) ===
// p' = R * (p + trans) (先平动后旋转)
tra_Rodvec(&p[0], coords, &rod_rig_F, trans_VecF)
// 实现: u[i] = Sum_j R[i*3+j]*(v[j] + t[j])
```

> **区别**: `Rodvec_tra` 用于目标自身姿态(先绕质心旋转再平移);
> `tra_Rodvec` 用于运动轨迹叠加(先平移再绕指定轴旋转)。

---

## 四、核心: 静态单站 V-ACA 算法 (`ACA_V.cpp`)

### 4.1 函数签名

```cpp
namespace ACA_v {
void ACA_v_calculate(
    int maxele,                    // 单个目标 RWG 边数
    int nums,                      // 总边数 = numTargets * maxele
    double th0, double dth, int nth,  // theta 参数
    double ph0, double dph, int nph,  // phi 参数
    vector<COMPLEX>& U_aca,        // <-- 输出 U 矩阵 (空间分量)
    vector<COMPLEX>& V_aca,        // <-- 输出 V 矩阵 (角度分量)
    Mesh& mesh, vector<Trianinfo>& Triangles, Points& points,
    COMPLEX k, COMPLEX eta, Edge_info Edges,
    vector<double*> trans_VecF,    // 各目标平动
    vector<double*> rot_VecF,      // 各目标旋转轴
    vector<double> rot_Angle,      // 各目标旋转角
    double pol_i);
}
```

### 4.2 算法伪代码

```
输入: 激励矩阵 B(edg1s * edg2s), edg1s=总边数, edg2s=角度数
输出: U(edg1s * rank), V(edg2s * rank) 使得 B ~ U * V^T

1. 初始化: index=0, I_k=0(首行), row_index={0}

2. 第一轮: 固定行 I_k=0, 遍历所有列
   for j in 0..edg2s-1:
       计算 B(0, j) -> v_in[j]
   选最大 |v_in[j]| -> J_k
   V 归一化: V(0,:) = v_in / v_in[J_k]

3. 第一轮: 固定列 J_k, 遍历所有行
   for i in 0..edg1s-1:
       计算 B(i, J_k) -> u_in[i]
   U(:,0) = u_in
   选最大 |u_in[i]| -> I_k (下一行)
   index = 1

4. 迭代 (while index < min(edg1s,edg2s)):
   a) 选新列: 固定行 I_k, 计算所有列
      v_in[j] = B(I_k,j) - Sum_{p=0}^{index-1} U(I_k,p)*V(j,p)
      选最大 |v_in[j]| -> J_k
      V(index,:) = v_in / v_in[J_k]

   b) 选新行: 固定列 J_k, 计算所有行
      u_in[i] = B(i,J_k) - Sum_{p=0}^{index-1} U(i,p)*V(J_k,p)
      U(:,index) = u_in
      选最大 |u_in[i]| -> I_k

   c) 更新 Frobenius 范数 ||B_approx||^2_F
      收敛: ||B_approx||^2 * eps^2 >= ||u_new||^2 * ||v_new||^2
      若满足 -> break

输出 rank = index
```

### 4.3 收敛判据

```cpp
if (Z_simlar_pingf * pow(eps_ACA, 2) >= U_f * V_f) break;
```

- Z_simlar_pingf = 已逼近的 Frobenius 范数平方 ||B_k||^2_F
- U_f * V_f = 新秩-1 分量的贡献 ||u_new||^2 * ||v_new||^2
- eps_ACA = 10^-3: 当新贡献的相对增量 < 0.1% 时停止

### 4.4 完整代码(逐行中文注释)

```cpp
namespace ACA_v {
void ACA_v_calculate(int maxele, int nums, ...) {
    th0 += 180.0;                          // 单站: 入射角 + 180 -> 散射角
    int I_k = 0, J_k, index = 0;           // 当前行/列/已完成秩数
    double Z_simlar_pingf = 0.0;           // ||B_approx||^2_F
    vector<int> row_index, col_index;      // 已选行/列
    int edg1s = nums;                      // 行数 = 总边数
    int edg2s = (nth==1) ? nph : nth;      // 列数 = 角度数
    double U_f, V_f, UU_VV_k, UU_k, VV_k;

    vector<COMPLEX> u_in, v_in;
    vector<COMPLEX> U_V_k(edg1s), U_k_V(edg2s);
    row_index.push_back(0);                // 首行 = 第0条边

    // ======== 第一轮: 选列 (固定行0, 找最大列) ========
    COMPLEX max_v(0), max_u(0);
    for (int s = 0; s < edg2s; ++s) {
        V_E[0]=V_H[0]=0;
        int k_tar = I_k / maxele;          // 该边属于哪个目标
        double angle = (nth==1) ? ph0+s*dph : th0+s*dth;

        vector<double> rod_rig_F(9);
        Rodrig(&rod_rig_F[0], rot_VecF[k_tar], rot_Angle[k_tar]);

        // 计算该边在该角度下的激励
        assemble_exic_vector(V_E, V_H, mesh, Triangles, points,
            I_k, k, eta0, ..., trans_VecF[k_tar], rod_rig_F, maxele);

        // CFIE 组合
        v_in.push_back(Aleph*V_E[0] + (1.0-Aleph)*V_H[0]);
    }

    // 选最大值列
    for (int i = 0; i < edg2s; i++)
        if (b_s(max_v, v_in[i])) { max_v = v_in[i]; J_k = i; }
    col_index.push_back(J_k);

    // 归一化 V 的第一行
    for (int i = 0; i < edg2s; i++) {
        V_aca.push_back(v_in[i] / max_v);
        V_f += F_n2(V_aca[i]);             // F_n2(x) = |x|^2
    }

    // ======== 第一轮: 选行 (固定列 J_k, 找最大行) ========
    for (int f = 0; f < edg1s; ++f) {
        V_E[0]=V_H[0]=0;
        int k_tar = f / maxele;
        double J_angle = (nth==1) ? ph0+J_k*dph : th0+J_k*dth;
        vector<double> rod_rig_F(9);
        Rodrig(&rod_rig_F[0], rot_VecF[k_tar], rot_Angle[k_tar]);

        assemble_exic_vector(V_E, V_H, mesh, Triangles, points,
            f, k, eta0, ..., trans_VecF[k_tar], rod_rig_F, maxele);

        u_in.push_back(Aleph*V_E[0] + (1.0-Aleph)*V_H[0]);
        U_aca.push_back(u_in[f]);
        U_f += F_n2(U_aca[f]);
    }

    Z_simlar_pingf += U_f * V_f;           // ||B||^2 += ||u1||^2 * ||v1||^2

    // 选下一行(未选过的)
    for (int i = 0; i < edg1s; i++)
        if (b_s(max_u, u_in[i]) && not_selected(row_index, i))
            { max_u = u_in[i]; I_k = i; }
    index++; row_index.push_back(I_k);

    // ======== 迭代主循环 ========
    while (index < min(edg1s, edg2s)) {
        U_f = V_f = 0;
        fill(u_in,v_in,U_V_k,U_k_V, 0);

        // --- 选新列 ---
        for (int s = 0; s < edg2s; ++s) {
            // [同上计算激励...]
            v_in[s] = Aleph*V_E[0] + (1.0-Aleph)*V_H[0];

            // *** 关键: 减去之前秩的贡献 ***
            for (int p = 0; p < index; p++)
                U_k_V[s] += U_aca[p*edg1s + I_k] * V_aca[p*edg2s + s];
            v_in[s] -= U_k_V[s];

            if (b_s(max_v, v_in[s]) && not_selected(col_index, s))
                { max_v = v_in[s]; J_k = s; }
        }
        col_index.push_back(J_k);

        // 归一化新V行
        for (int i = 0; i < edg2s; i++) {
            V_aca.push_back(v_in[i] / max_v);
            V_f += F_n2(V_aca[index*edg2s + i]);
        }

        // --- 选新行 ---
        for (int f = 0; f < edg1s; ++f) {
            // [同上计算激励...]
            u_in[f] = Aleph*V_E[0] + (1.0-Aleph)*V_H[0];

            // *** 关键: 减去之前秩的贡献 ***
            for (int p = 0; p < index; p++)
                U_V_k[f] += U_aca[p*edg1s + f] * V_aca[p*edg2s + J_k];
            u_in[f] -= U_V_k[f];

            U_aca.push_back(u_in[f]);
            U_f += F_n2(U_aca[index*edg1s + f]);

            if (b_s(max_u, u_in[f]) && not_selected(row_index, f))
                { max_u = u_in[f]; I_k = f; }
        }

        // --- 更新 Frobenius 范数 ---
        for (int p = 0; p < index; p++) {
            UU_k = VV_k = 0;
            for (int i = 0; i < edg1s; i++)
                UU_k += F_n1(U_aca[p*edg1s+i] * U_aca[index*edg1s+i]);
            for (int i = 0; i < edg2s; i++)
                VV_k += F_n1(V_aca[index*edg2s+i] * V_aca[p*edg2s+i]);
            UU_VV_k += UU_k * VV_k;
        }
        Z_simlar_pingf += U_f * V_f + 2.0 * UU_VV_k;

        index++;
        // --- 收敛判断 ---
        if (Z_simlar_pingf * pow(eps_ACA, 2) >= U_f * V_f) break;
        row_index.push_back(I_k);
    }
    cout << "angles = " << index << endl;
}
}
```

---

## 五、核心: 运动目标 V-ACA 算法 (`ACA_v_mov.cpp`)

### 5.1 与静态版本的差异

| 项目 | 静态 ACA_V | 运动 ACA_v_mov |
|------|-----------|---------------|
| 列含义 | 入射角度 | 时间步 |
| 列数 `edg2s` | nth 或 nph | time_mov_s |
| 入射方向 | 变角度扫掠 | 固定 (th0, ph0) |
| 坐标变换 | Rodrigues(目标姿态) | 双轴 Rodrigues + 平动 |
| 旋转变换数 | 1 个矩阵 | 3 个矩阵相乘 |

### 5.2 运动变换构建(每时间步)

```cpp
// 1. 平动: dx = t * vx,  dy = t * vy,  dz = t * vz
trans_F[0] = s_edg * trans_t_s[0];
trans_F[1] = s_edg * trans_t_s[1];
trans_F[2] = s_edg * trans_t_s[2];

// 2. 旋转角: 匀速转动
ro_angle_t1 = s_edg * rot_ang_s1;   // 绕轴1的累积转角
ro_angle_t2 = s_edg * rot_ang_s2;   // 绕轴2的累积转角

// 3. R1 = Rodrigues(axis1, angle1)
Rodrig(&rod_rig_1[0], &rot_t_s1[0], ro_angle_t1);

// 4. axis2' = R1 * axis2 (第二轴受第一旋转影响)
Rodvec(&rot_t_s21[0], &rot_t_s2[0], &rod_rig_1[0]);

// 5. R2 = Rodrigues(axis2', angle2)
Rodrig(&rod_rig_2[0], &rot_t_s21[0], ro_angle_t2);

// 6. R_total = R2 * R1
for (int i=0; i<3; i++)
    for (int j=0; j<3; j++)
        for (int k=0; k<3; k++)
            rod_rig_3[i*3+j] += rod_rig_2[i*3+k] * rod_rig_1[k*3+j];

// 7. 联合变换计算激励
assemble_exic_vector(V_E, V_H, ..., I_k, k, eta0, th0, ph0, pol_i, Edges,
                     &trans_F[0], rod_rig_3, maxele);
```

### 5.3 完整代码

```cpp
namespace ACA_v_mov {
void ACA_v_mov_cal(int maxele, int nums, double th0, double ph0,
    vector<COMPLEX>& U_aca, vector<COMPLEX>& V_aca,
    Mesh& mesh, vector<Trianinfo>& Triangles, Points& points,
    COMPLEX k, COMPLEX eta, Edge_info Edges,
    vector<double> trans_t_s, vector<double> rot_t_s1,
    double rot_ang_s1, int time_mov_s, double pol_i,
    vector<double> rot_t_s2, double rot_ang_s2)
{
    // ======== 初始化 ========
    th0 += 180.0;
    int I_k = 0, J_k, index = 0;
    double Z_simlar_pingf = 0.0;
    vector<int> row_index, col_index;
    int edg1s = nums;              // 行 = 总边数
    int edg2s = time_mov_s;        // *** 列 = 时间步数 ***

    vector<COMPLEX> u_in, v_in, U_V_k(edg1s), U_k_V(edg2s);
    row_index.push_back(0);

    // ======== 第一轮: 选列 ========
    COMPLEX max_v(0), max_u(0);
    for (int s = 0; s < edg2s; ++s) {
        // *** 构建时间步 s 的变换 ***
        vector<double> trans_F(3);
        trans_F[0]=s*trans_t_s[0]; trans_F[1]=s*trans_t_s[1];
        trans_F[2]=s*trans_t_s[2];
        double a1=s*rot_ang_s1, a2=s*rot_ang_s2;
        vector<double> R1(9), R2(9), Rtotal(9), axis2r(3);
        Rodrig(&R1[0], &rot_t_s1[0], a1);
        Rodvec(&axis2r[0], &rot_t_s2[0], &R1[0]);
        Rodrig(&R2[0], &axis2r[0], a2);
        // Rtotal = R2 * R1
        for (int i=0;i<3;i++)
            for (int j=0;j<3;j++)
                for (int k=0;k<3;k++)
                    Rtotal[i*3+j] += R2[i*3+k]*R1[k*3+j];

        // 计算固定入射方向的激励
        V_E[0]=V_H[0]=0;
        assemble_exic_vector(V_E,V_H,mesh,Triangles,points,
            I_k,k,eta0,th0,ph0,pol_i,Edges,&trans_F[0],Rtotal,maxele);
        v_in.push_back(Aleph*V_E[0] + (1.0-Aleph)*V_H[0]);
    }

    // 选最大列, 归一化 V
    for (int i=0;i<edg2s;i++)
        if (b_s(max_v,v_in[i])) { max_v=v_in[i]; J_k=i; }
    col_index.push_back(J_k);
    for (int i=0;i<edg2s;i++) {
        V_aca.push_back(v_in[i]/max_v);
        V_f += F_n2(V_aca[i]);
    }

    // ======== 第一轮: 选行 (固定时间步 J_k) ========
    for (int f=0; f<edg1s; ++f) {
        // [构建 J_k 时刻的变换, 同上]
        V_E[0]=V_H[0]=0;
        assemble_exic_vector(V_E,V_H,mesh,Triangles,points,
            f,k,eta0,th0,ph0,pol_i,Edges,&trans_F[0],Rtotal,maxele);
        u_in.push_back(Aleph*V_E[0] + (1.0-Aleph)*V_H[0]);
        U_aca.push_back(u_in[f]);
        U_f += F_n2(U_aca[f]);
    }
    Z_simlar_pingf += U_f * V_f;
    for (int i=0;i<edg1s;i++)
        if (b_s(max_u,u_in[i]) && not_selected(row_index,i))
            { max_u=u_in[i]; I_k=i; }
    index++; row_index.push_back(I_k);

    // ======== 迭代主循环 ========
    while (index < min(edg1s, edg2s)) {
        // [结构同静态版本, 区别仅在于列遍历用 time_mov_s
        //  且每次计算激励前构建完整运动变换]
        // ... 收敛条件相同 ...
    }
    cout << "angles = " << index << endl;
}
}
```

### 5.4 多目标扩展 `ACA_v_mov_cal_bua`

在有多个目标时 `_bua` 版本叠加两层变换:

```cpp
// 目标 k 在时间步 s 的完整位置
trans_F[0] = s*trans_t_s[0] + trans_VecF[k][0];   // 轨迹 + 目标偏移
// 旋转: rod_rig_F(目标自身姿态) * R_total(运动)
// rod_rig_4 = rod_rig_3 * rod_rig_2  (多一层矩阵乘)
```

---

## 六、主程序集成: 完整调用流程

文件 `PMCHWT.cpp`，条件: `Moving==1 && num_Targets==1 && Dan_aca_V==1` (L905-L970)。

### 6.1 步骤 1: ACA 分解右边向量矩阵

```cpp
Btree nodeV = new Bnode;
ACA_v_mov::ACA_v_mov_cal(
    maxele,                              // 每目标边数
    num_Targets * maxele,                // 总边数 = 行数
    th0, ph0,                            // 固定入射方向
    nodeV->U_matrix,                     // <-- U(edges * rank)
    nodeV->V_matrix,                     // <-- V(timesteps * rank)
    mesh, Triangles, points, k0, eta0, Edges,
    trans_t_s, rot_t_s1, rot_ang_s1,     // 运动参数
    time_mov_s, pol_i,
    rot_t_s2, rot_ang_s2);
```

### 6.2 步骤 2: 求解 rank 个等效线性系统

```cpp
int rank = nodeV->U_matrix.size() / (num_Targets * maxele);
vector<COMPLEX> U_new;                    // 等效解矩阵

for (int m = 0; m < rank; m++) {
    // 以 U 的第 m 列作为右边向量
    COMPLEX* x_est = new COMPLEX[num_Targets * maxele]();
    COMPLEX* e_R   = new COMPLEX[num_Targets * maxele];
    for (int i = 0; i < num_Targets * maxele; i++)
        e_R[i] = nodeV->U_matrix[m * num_Targets * maxele + i];

    // MGMRES 迭代: Z * x_m = u_m
    mgmres_stnew(num_Targets * maxele,
        num_Targets*num_Targets*maxele*maxele,
        x_est, e_R, itr_max, mr, tol_abs, tol_rel,
        near_col, near_row, value_near,
        U_V_p, Cude_index, cude_Distan, num_Targets,
        index_X, T_rcv, T_edg_near, rtt);

    // 保存等效解
    for (int i = 0; i < num_Targets * maxele; i++)
        U_new.push_back(x_est[i]);
    delete[] x_est; delete[] e_R;
}
```

### 6.3 步骤 3: 合成每个时间步的解并计算 RCS

```cpp
th = th0; ph = ph0;
for (int t_s = 0; t_s < time_mov_s; t_s++) {
    // --- 线性叠加: x(t_s) = Sum_m x_m * V(m, t_s) ---
    COMPLEX* x_sol = new COMPLEX[num_Targets * maxele]();
    for (int i = 0; i < num_Targets * maxele; i++)
        for (int m = 0; m < rank; m++)
            x_sol[i] += U_new[m * N + i]
                      * nodeV->V_matrix[m * time_mov_s + t_s];

    // --- 构建 t_s 时刻变换 ---
    vector<double> trans_F(3);
    trans_F[0]=t_s*trans_t_s[0]; trans_F[1]=t_s*trans_t_s[1];
    trans_F[2]=t_s*trans_t_s[2];
    double a1=t_s*rot_ang_s1, a2=t_s*rot_ang_s2;
    vector<double> R1(9),R2(9),Rtotal(9),axis2r(3);
    Rodrig(&R1[0],&rot_t_s1[0],a1);
    Rodvec(&axis2r[0],&rot_t_s2[0],&R1[0]);
    Rodrig(&R2[0],&axis2r[0],a2);
    for(int i=0;i<3;i++)for(int j=0;j<3;j++)for(int k=0;k<3;k++)
        Rtotal[i*3+j]+=R2[i*3+k]*R1[k*3+j];

    // --- 计算 RCS ---
    vector<COMPLEX> RCS_H(3), RCS_E(3);
    for (int num = 0; num < num_Targets; num++)
        RCS::calculate_rcs(x_sol, mesh, Triangles, points,
            maxele, k0, eta0, th, ph, RCS_H, RCS_E, pol_i, Edges,
            num, &trans_F[0], Rtotal);

    double rcs = 10.0 * log10( (1.0/(4.0*PI))
        * pow(k0.real()*eta0.real(),2)
        * (F_n2(RCS_H[0]+RCS_E[0])
         + F_n2(RCS_H[1]+RCS_E[1])
         + F_n2(RCS_H[2]+RCS_E[2])) );
    ofs << th << " " << ph << " " << fre << " " << rcs << endl;
    delete[] x_sol;
}
```

### 6.4 线性叠加原理

```
已知: Z * x = b

ACA 将 b 矩阵分解:  B ~ U * V^T
即:  B(:,t) ~ Sum_m U(:,m) * V(t,m)

对每个模态 m, 求解:  Z * x_m = U(:,m)     (MGMRES)

则任意时刻 t 的解:
  x(t) = Sum_m x_m * V(t,m)               (纯代数, 极快)

=> 只需 rank(~15) 次求解, 而非 time_mov_s(~2000) 次!
```

### 6.5 加速效果

| 参数 | 典型值 |
|------|--------|
| 时间步数 | 2000 |
| ACA 秩 | ~15 |
| 求解加速比 | 2000/15 ~ 133x |
| 额外开销 | rank 次 MGMRES + 时间步合成(可忽略) |

---

## 七、Rodrigues 旋转公式 (`Products.h`)

```cpp
// 公式: R = cos(theta)*I + (1-cos(theta))*u*u^T + sin(theta)*[u]x
// 输入: uvec(旋转轴单位向量, 3元素), CONST(旋转角, 度)
// 输出: res(3*3 旋转矩阵, 按行展开为 9 元素)
inline void Rodrig(T* res, T* uvec, T& CONST) {
    double c = cos(CONST * PI / 180.0);
    double s = sin(CONST * PI / 180.0);
    double t = 1.0 - c;

    res[0] = c + uvec[0]*uvec[0]*t;         // R11
    res[1] = uvec[0]*uvec[1]*t - uvec[2]*s; // R12
    res[2] = uvec[0]*uvec[2]*t + uvec[1]*s; // R13
    res[3] = uvec[1]*uvec[0]*t + uvec[2]*s; // R21
    res[4] = c + uvec[1]*uvec[1]*t;         // R22
    res[5] = uvec[1]*uvec[2]*t - uvec[0]*s; // R23
    res[6] = uvec[2]*uvec[0]*t - uvec[1]*s; // R31
    res[7] = uvec[2]*uvec[1]*t + uvec[0]*s; // R32
    res[8] = c + uvec[2]*uvec[2]*t;         // R33
}

// 向量乘法
inline void Rodvec(T* u, T* v, T* R) {
    // u = R * v (仅旋转)
    for (int i=0; i<3; i++)
        for (int j=0; j<3; j++)
            u[i] += R[i*3+j] * v[j];
}
inline void Rodvec_tra(T* u, T* v, T* R, T* t) {
    // u = R * v + t (旋转+平动)
    for (int i=0; i<3; i++) {
        for (int j=0; j<3; j++) u[i] += R[i*3+j] * v[j];
        u[i] += t[i];
    }
}
inline void tra_Rodvec(T* u, T* v, T* R, T* t) {
    // u = R * (v + t) (平动+旋转)
    for (int i=0; i<3; i++)
        for (int j=0; j<3; j++)
            u[i] += R[i*3+j] * (v[j] + t[j]);
}
```

辅助工具函数:

| 函数 | 功能 | 公式 |
|------|------|------|
| `b_s(a, b)` | 复数模比较 | `|a|^2 < |b|^2` |
| `F_n2(x)` | 模平方 | `Re(x)^2 + Im(x)^2` |
| `F_n1(x)` | 模 | `sqrt(Re(x)^2 + Im(x)^2)` |
| `dot(u, v)` | 向量内积 | `u1*v1 + u2*v2 + u3*v3` |
| `cross(r, u, v)` | 向量叉乘 | `r = u * v` |

---

## 八、总结

本项目将 ACA 低秩分解从传统的阻抗矩阵 Z 扩展到 **右边激励向量矩阵**，核心思想:

1. **构建激励矩阵**: 边缘索引为行, 角度/时间步为列
2. **ACA 贪婪分解**: 交替选择最大残差的行和列, 逐步低秩逼近
3. **残差更新**: 每轮计算 = 精确值 - 之前秩的贡献之和
4. **线性叠加**: 利用 MoM 的线性性, rank 次求解替代 thousands 次
5. **运动延拓**: 通过 Rodrigues 双轴旋转 + 匀速平动, 推广到运动目标

**完整数据流**:

```
[网格数据 + 运动参数(trans_t_s, rot_t_s1/2, rot_ang_s1/2)]
    |
    v
[ACA 分解激励矩阵] --> U(edges * rank) + V(timesteps * rank)
    |                        |
    v                        v
[Z * x_m = U(:,m)]    [保存 V 矩阵]
    |   (MGMRES * rank)
    v
[U_new: 等效解矩阵]
    |
    v
[每个时间步: x(t) = Sum_m U_new(:,m) * V(m,t)]
    |
    v
[计算 RCS * 输出]
```

**关键文件索引**:

| 文件 | 行号/内容 |
|------|----------|
| `PMCHWT.cpp` | L918 运动ACA-V入口, L905-970 完整流程 |
| `PMCHWT.cpp` | L690 静态ACA-V入口, L677-745 完整流程 |
| `ACA_v_mov.cpp` | 全文: 运动目标右边向量ACA分解 |
| `ACA_V.cpp` | 全文: 静态目标右边向量ACA分解 |
| `ACA_v_mov.h` | 声明 `ACA_v_mov_cal` 和 `ACA_v_mov_cal_bua` |
| `excvecE.cpp` | L14 `assemble_exic_vector`, L105 `assemble_exic_vector_mov` |
| `mymathparater.h` | Aleph=0.2, eps_ACA=1e-3, eta0=377 Ohm |
| `Products.h` | Rodrigues旋转, 向量运算(F_n2, b_s, cross, dot等) |
| `Z_MoM.cpp` | PMCHWT组合: Aleph*Z_EFIE + (1-Aleph)*Z_MFIE |
| `EMsource.cpp/h` | 平面波源: initPW(极化, theta, phi) |

---

> 文档生成: 2026-06-03
> 项目路径: `D:\CEM\MLH\ACA-MOM-finally\pmchw_fortri`
