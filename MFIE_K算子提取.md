# MFIE 中 K 算子的公式与代码映射及完整提取

## 1. 范围

本文仅分析并提取 `Code/MFIE.h` 与 `Code/MFIE.cpp` 中的 MFIE/K 算子代码，公式依据为 `任务1.md`，不展开其他代码文件。

## 2. 参考公式

MFIE 为

```math
\frac{1}{2}\bar{\mathbf J}_s+
\hat{\mathbf n}\times\mathcal K(\bar{\mathbf J}_s)
=
\hat{\mathbf n}\times\mathbf H^i(\mathbf r).
```

K 算子及格林函数为

```math
\mathcal K(\mathbf X)
=
\iint_{S_0}\mathbf X(\mathbf r')\times
\nabla G_0(\mathbf r,\mathbf r')\,\mathrm dS',
\qquad
G_0=\frac{e^{-\mathrm j k_0R}}{4\pi R},
\quad R=|\mathbf r-\mathbf r'|.
```

令 `Rvec = r-r'`，则

```math
\nabla_{\mathbf r}G_0
=
-\frac{\mathbf R}{4\pi R^3}(1+\mathrm jk_0R)e^{-\mathrm jk_0R}.
```

利用叉乘反对称性，

```math
\mathbf f_n\times\nabla G_0
=
\frac{(1+\mathrm jk_0R)e^{-\mathrm jk_0R}}{4\pi R^3}
(\mathbf R\times\mathbf f_n).
```

这正是代码中 `gradG`、`Rvec` 与第一次叉乘的来源。

## 3. RWG 伽辽金离散

在每个 RWG 支撑三角形上，

```math
\mathbf f(\mathbf r)
=
s\frac{l}{2A}(\mathbf r-\mathbf r_0),
\qquad s\in\{-1,+1\}.
```

K 算子部分的测试矩阵元为

```math
K_{mn}
=
\iint
\mathbf f_m(\mathbf r)\cdot
\left[
\hat{\mathbf n}(\mathbf r)\times
\left(
\mathbf f_n(\mathbf r')\times\nabla G_0
\right)
\right]
\,\mathrm dS'\mathrm dS.
```

代码用 `pomocvF=r-r_m^0`、`pomocvS=r'-r_n^0` 保存未除以 `2A` 的局部 RWG 向量。三角形高斯积分的面积因子与两个 `1/(2A)` 相消；边长、RWG 正负号和格林函数的 `1/(4π)` 最后由

```cpp
(lF * lS) / (4.0 * Pi) * signSymbolF * signSymbolS
```

统一乘入。

## 4. 公式与代码映射

| 数学对象/步骤 | 代码实现 |
|---|---|
| MFIE 矩阵元入口 | `MFIE::computeMFIE_Zij` |
| K 算子奇异补偿 | `MFIE::singularityMFIE` |
| 场/源 RWG 的正负支撑三角形 | `fi=0..1`、`sj=0..1` |
| 局部 RWG 系数 `1/(2A)` | `rec2Af`、`rec2As` |
| 场点、源点及 `R=r-r'` | `fieldPoint`、`sourcePoint`、`Rvec` |
| 去掉 `1/(4π)` 的梯度核系数 | `gradG=exp(-J*k_*R)*(1+J*k_*R)/R^3` |
| `R×f_n` | `cross(CROSS1, Rvec, pomocvS)` |
| `n×(R×f_n)` | `cross(CROSS2, normalVectorF, CROSS1)` |
| 以场 RWG 作伽辽金测试 | `dot(pomocvF, CROSS2)` |
| MFIE 的 `1/2` 跳跃项 | 同面分支中的单重 `numI` 积分 |
| RWG 边长及正负号 | `lF/lS`、`signSymbolF/signSymbolS` |
| `1/(4π)` | 最终累加 `Zij_near_M` 时除以 `4*Pi` |
| 同面核的平滑余项 | 同面分支中的复数变量 `G` |
| `1/R^3` 与 `k^2/(2R)` 解析项 | `singularityMFIE` 中的 `Ivec1R3/Isca1R3` 与 `Ivec1R/Isca1R` |

## 5. 同面奇异核与跳跃项

梯度核的标量部分在 `R→0` 时展开为

```math
\frac{e^{-\mathrm jkR}(1+\mathrm jkR)}{R^3}
=
\frac{1}{R^3}
+
\frac{k^2}{2R}
+
\text{有限余项}.
```

所以同面分支使用

```cpp
G = exp(-J*k_*R)*(1+J*k_*R)/R^3
    - 1/R^3
    - 0.5*k_*k_/R;
```

对有限余项作双重高斯积分；`singularityMFIE` 分别解析计算被减去的 `1/R^3` 和 `k^2/(2R)` 部分，再通过 `alok` 加回。

MFIE 方程还包含独立的 `1/2` 跳跃项。代码在同面分支以单重高斯积分计算

```math
\frac12\int_{S_m}\mathbf f_m\cdot\mathbf f_n\,\mathrm dS.
```

局部先乘 `4π`，在最终统一除以 `4π` 后得到正确的跳跃项矩阵元。

## 6. 提取后的完整头文件

```cpp
// MFIE.h
#pragma once
#ifndef MFIE_H
#define MFIE_H

#include "MyStruct.h"
#include "GaussPoints.h"

namespace MFIE {
	void computeMFIE_Zij(const RWGBase* rwgField, const RWGBase* rwgSource, const std::complex<double> k_,
		const gaussPoints& gp, std::complex<double>& Zij_near_M);

	void singularityMFIE(double& areaTriF, double& areaTriS,
		double* vertexF1, double* vertexF2, double* vertexF3,
		double* nonPublicEdge_vertexF, double* normalVectorF, const double recAf,
		double* vertexS1, double* vertexS2, double* vertexS3,
		double* nonPublicEdge_vertexS, double* normalVectorS, const double recAs,
		const gaussPoints& gp, const std::complex<double> k_, std::complex<double>& alok);
}
#endif // !MFIE_H
```

## 7. 提取后的完整实现文件

`computeMFIE_Zij` 负责 K 算子常规装配、同面核拆分和 `1/2` 跳跃项；`singularityMFIE` 负责 `1/R^3` 与 `1/R` 奇异部分的解析补偿。

```cpp
// MFIE.cpp
#include "MFIE.h"
#include "OverloadAlgo.h"

namespace MFIE {
	void computeMFIE_Zij(const RWGBase* rwgField, const RWGBase* rwgSource, const std::complex<double> k_,
		const gaussPoints& gp, std::complex<double>& Zij_near_M)
	{
		int rwgFieldTriId, rwgSourceTriId;
		double areaTriF, areaTriS;
		double rec2Af, rec2As;
		double* nonPublicEdge_vertexF, * nonPublicEdge_vertexS;
		double* vertexF1, * vertexF2, * vertexF3;
		double* vertexS1, * vertexS2, * vertexS3;
		double* normalVectorF;
		double* normalVectorS;
		double lF, lS;

		std::vector<double> fieldPoint(3), pomocvF(3), RWGf(3);
		std::vector<double> sourcePoint(3), pomocvS(3), RWGs(3);

		std::vector<double> CROSS1(3), CROSS2(3);

		std::vector<double> Rvec(3);
		std::vector<double> Rivec(3);

		for (int fi = 0; fi < 2; ++fi) {
			lF = rwgField->edgeLength;
			if (fi == 0) {
				rwgFieldTriId = rwgField->negativeFace->faceId;
				// ��RWG�����������������
				areaTriF = rwgField->negativeFace->triangularArea;
				rec2Af = 1.0 / (2.0 * areaTriF);
				// ע�⣺����ķǹ����߶����Ǹ������εĶ���
				nonPublicEdge_vertexF = &(rwgField->freeVertexNegative->x);
				// ��ó�RWG������������������
				vertexF1 = &(rwgField->negativeFace->vertex1->x);
				vertexF2 = &(rwgField->negativeFace->vertex2->x);
				vertexF3 = &(rwgField->negativeFace->vertex3->x);
				normalVectorF = &(rwgField->negativeFace->externalNormalVector.x);
			}
			else {
				rwgFieldTriId = rwgField->positiveFace->faceId;
				// ��RWG�����������������
				areaTriF = rwgField->positiveFace->triangularArea;
				rec2Af = 1.0 / (2.0 * areaTriF);
				// ��ȡ��RWG�������ķǹ����߶���
				nonPublicEdge_vertexF = &(rwgField->freeVertexPositive->x);
				// ��ó�RWG������������������
				vertexF1 = &(rwgField->positiveFace->vertex1->x);
				vertexF2 = &(rwgField->positiveFace->vertex2->x);
				vertexF3 = &(rwgField->positiveFace->vertex3->x);
				normalVectorF = &(rwgField->positiveFace->externalNormalVector.x);
			}

			for (int sj = 0; sj < 2; ++sj) {
				lS = rwgSource->edgeLength;
				if (sj == 0) {
					rwgSourceTriId = rwgSource->negativeFace->faceId;
					// ԴRWG�����������������
					areaTriS = rwgSource->negativeFace->triangularArea;
					rec2As = 1.0 / (2.0 * areaTriS);
					// ��ȡ��RWG�������ķǹ����߶���
					// ע�⣺����ķǹ����߶����Ǹ������εĶ���
					nonPublicEdge_vertexS = &(rwgSource->freeVertexNegative->x);
					// ��ó�RWG������������������
					vertexS1 = &(rwgSource->negativeFace->vertex1->x);
					vertexS2 = &(rwgSource->negativeFace->vertex2->x);
					vertexS3 = &(rwgSource->negativeFace->vertex3->x);
					// ���ԴRWG�������������εķ�����
					normalVectorS = &(rwgSource->negativeFace->externalNormalVector.x);
				}
				else {
					rwgSourceTriId = rwgSource->positiveFace->faceId;
					//  ԴRWG�����������������
					areaTriS = rwgSource->positiveFace->triangularArea;
					rec2As = 1.0 / (2.0 * areaTriS);
					// ��ȡԴRWG�������ķǹ����߶���
					nonPublicEdge_vertexS = &(rwgSource->freeVertexPositive->x);
					// ���ԴRWG������������������
					vertexS1 = &(rwgSource->positiveFace->vertex1->x);
					vertexS2 = &(rwgSource->positiveFace->vertex2->x);
					vertexS3 = &(rwgSource->positiveFace->vertex3->x);
					// ���ԴRWG�������������εķ�����
					normalVectorS = &(rwgSource->positiveFace->externalNormalVector.x);
				}

				std::complex<double> alok(0.0, 0.0);
				std::complex<double> alok1(0.0, 0.0);

				if (rwgFieldTriId != rwgSourceTriId) {
					// ���ù�ʽ: $\sum_{f=1}^{n_P}\sum_{s=1}^{n_P}w_fw_s\mathbf{f}_m(\mathbf{r}_f)\cdot\left[\hat{n}_f\times((\mathbf{r}_f-\mathbf{r}_s)/R^3(1+ikR)e^{-ikR}\times\mathbf{f}_n(\mathbf{r}_s))\right]$
					// ѭ����������������Ԫ�ϵĻ��ֵ�
					for (int numI = 0; numI < gp.N_points_; ++numI) {
						// ���㳡RWG��������������Ԫ�ϵĻ��ֵ�����
						fieldPoint[0] = gp.l1[numI] * vertexF1[0] + gp.l2[numI] * vertexF2[0] + gp.l3[numI] * vertexF3[0];
						fieldPoint[1] = gp.l1[numI] * vertexF1[1] + gp.l2[numI] * vertexF2[1] + gp.l3[numI] * vertexF3[1];
						fieldPoint[2] = gp.l1[numI] * vertexF1[2] + gp.l2[numI] * vertexF2[2] + gp.l3[numI] * vertexF3[2];
						// ��ȡ�ó����ֵ�Ȩ��
						double wF = gp.weight[numI];

						// ѭ������Դ��������Ԫ�ϵĻ��ֵ�
						for (int numJ = 0; numJ < gp.N_points_; ++numJ) {
							// ����ԴRWG��������������Ԫ�ϵĻ��ֵ�����
							sourcePoint[0] = gp.l1[numJ] * vertexS1[0] + gp.l2[numJ] * vertexS2[0] + gp.l3[numJ] * vertexS3[0];
							sourcePoint[1] = gp.l1[numJ] * vertexS1[1] + gp.l2[numJ] * vertexS2[1] + gp.l3[numJ] * vertexS3[1];
							sourcePoint[2] = gp.l1[numJ] * vertexS1[2] + gp.l2[numJ] * vertexS2[2] + gp.l3[numJ] * vertexS3[2];
							// ��ȡ��Դ���ֵ�Ȩ��
							double wS = gp.weight[numJ];
							// ����ѡ����RWG��������������Ԫ�ϵĻ��ֵ㵽ԴRWG��������������Ԫ�ϵĻ��ֵ�ľ���
							// ��˹���ֹ���߻�ѡ��ͬһ��
							double R = norm(&fieldPoint[0], &sourcePoint[0]);

							sub(&Rvec[0], &fieldPoint[0], &sourcePoint[0]);

							std::complex<double> gradG = (exp(-J * k_ * R) / (R * R * R)) * (1.0 + J * k_ * R);

							sub(&pomocvF[0], &fieldPoint[0], nonPublicEdge_vertexF);
							sub(&pomocvS[0], &sourcePoint[0], nonPublicEdge_vertexS);

							cross(&CROSS1[0], &Rvec[0], &pomocvS[0]);
							cross(&CROSS2[0], normalVectorF, &CROSS1[0]);

							std::complex<double> contrib = wF * wS * (gradG * dot(&pomocvF[0], &CROSS2[0]));
							alok1 += contrib;
						}
					}
				}
				else {
					// ���ù�ʽ: $\sum_{f}\sum_{f}4w_{s}w_{s}A_{f}A_{s}\left[\nabla G(\mathbf{r}_{f},\mathbf{r}_{s})\cdot\left((\mathbf{r}_{f}-\mathbf{r}_{s})\times\mathbf{f}_{n}(\mathbf{r}_{s})\cdot(\mathbf{f}_{m}(\mathbf{r}_{f})\times\hat{n}_{f})\right)\right]$
					for (int numI = 0; numI < gp.N_points_; ++numI) {
						// ���㳡RWG��������������Ԫ�ϵĻ��ֵ�����
						fieldPoint[0] = gp.l1[numI] * vertexF1[0] + gp.l2[numI] * vertexF2[0] + gp.l3[numI] * vertexF3[0];
						fieldPoint[1] = gp.l1[numI] * vertexF1[1] + gp.l2[numI] * vertexF2[1] + gp.l3[numI] * vertexF3[1];
						fieldPoint[2] = gp.l1[numI] * vertexF1[2] + gp.l2[numI] * vertexF2[2] + gp.l3[numI] * vertexF3[2];
						// ��ȡ�ó����ֵ�Ȩ��
						double wF = gp.weight[numI];

						// ѭ������Դ��������Ԫ�ϵĻ��ֵ�
						for (int numJ = 0; numJ < gp.N_points_; ++numJ) {
							// ����ԴRWG��������������Ԫ�ϵĻ��ֵ�����
							sourcePoint[0] = gp.l1[numJ] * vertexS1[0] + gp.l2[numJ] * vertexS2[0] + gp.l3[numJ] * vertexS3[0];
							sourcePoint[1] = gp.l1[numJ] * vertexS1[1] + gp.l2[numJ] * vertexS2[1] + gp.l3[numJ] * vertexS3[1];
							sourcePoint[2] = gp.l1[numJ] * vertexS1[2] + gp.l2[numJ] * vertexS2[2] + gp.l3[numJ] * vertexS3[2];
							// ��ȡ��Դ���ֵ�Ȩ��
							double wS = gp.weight[numJ];

							// ����ѡ����RWG��������������Ԫ�ϵĻ��ֵ㵽ԴRWG��������������Ԫ�ϵĻ��ֵ�ľ���
							double R = norm(&fieldPoint[0], &sourcePoint[0]);

							sub(&Rvec[0], &fieldPoint[0], &sourcePoint[0]);

							std::complex<double> G;

							if (R < eps0) {
								G = std::complex<double>(0.0, 0.0);
							}
							else {
								G = (exp(-J * k_ * R) / (R * R * R)) * (1.0 + J * k_ * R) - 1.0 / (R * R * R) - 0.5 * (k_ * k_) * (1.0 / R);
							}

							sub(&pomocvF[0], &fieldPoint[0], nonPublicEdge_vertexF);
							sub(&pomocvS[0], &sourcePoint[0], nonPublicEdge_vertexS);

							sub(&Rivec[0], &fieldPoint[0], nonPublicEdge_vertexS);

							cross(&CROSS1[0], &Rivec[0], &pomocvS[0]);
							cross(&CROSS2[0], &pomocvF[0], normalVectorF);

							std::complex<double> contrib = wF * wS * (G * dot(&CROSS1[0], &CROSS2[0]));
							alok1 += contrib;
							//temp2 += contrib;
						}
					}
					// \frac{1}{2}\int_{S_{m}}\mathbf{\Lambda}_{m}\mathbf{\Lambda}_{n}dS_{m}
					// $\frac{1}{2}\sum_fw_f2A_{f}\left(\mathbf{f}_m(\mathbf{r}_f)\cdot\mathbf{f}_n(\mathbf{r}_f)\right)4\pi$
					for (int numI = 0; numI < gp.N_points_; ++numI) {
						// ���㳡RWG��������������Ԫ�ϵĻ��ֵ�����
						fieldPoint[0] = gp.l1[numI] * vertexF1[0] + gp.l2[numI] * vertexF2[0] + gp.l3[numI] * vertexF3[0];
						fieldPoint[1] = gp.l1[numI] * vertexF1[1] + gp.l2[numI] * vertexF2[1] + gp.l3[numI] * vertexF3[1];
						fieldPoint[2] = gp.l1[numI] * vertexF1[2] + gp.l2[numI] * vertexF2[2] + gp.l3[numI] * vertexF3[2];

						double wF = gp.weight[numI];

						sub(&pomocvF[0], &fieldPoint[0], nonPublicEdge_vertexF);
						mul(&RWGf[0], &pomocvF[0], rec2Af);

						sub(&pomocvS[0], &fieldPoint[0], nonPublicEdge_vertexS);
						mul(&RWGs[0], &pomocvS[0], rec2As);

						std::complex<double> contrib = 4.0 * Pi * wF * areaTriF * dot(&RWGf[0], &RWGs[0]);
						alok += contrib;
						//temp3 += contrib;
					}

					singularityMFIE(areaTriF, areaTriS,
						vertexF1, vertexF2, vertexF3,
						nonPublicEdge_vertexF, normalVectorF, rec2Af,
						vertexS1, vertexS2, vertexS3,
						nonPublicEdge_vertexS, normalVectorS, rec2As,
						gp, k_, alok);

					alok1 += alok;
				}

				// ��Ԫ�������������������Σ���Ϊ1������Ϊ-1
				const double signSymbolF = (fi == 0) ? -1.0 : 1.0;
				const double signSymbolS = (sj == 0) ? -1.0 : 1.0;
				Zij_near_M += (lF * lS) / (4.0 * Pi) * signSymbolF * signSymbolS * alok1;
			}
		}
	}

	void singularityMFIE(double& areaTriF, double& areaTriS,
		double* vertexF1, double* vertexF2, double* vertexF3,
		double* nonPublicEdge_vertexF, double* normalVectorF, const double recAf,
		double* vertexS1, double* vertexS2, double* vertexS3,
		double* nonPublicEdge_vertexS, double* normalVectorS, const double recAs,
		const gaussPoints& gp, const std::complex<double> k_, std::complex<double>& alok)
	{
		std::vector<double> fieldPoint(3);

		std::vector<double> npom1(3);
		std::vector<double> npom2(3);
		std::vector<double> npom3(3);
		std::vector<double> rho1(3);
		std::vector<double> rho2(3);
		std::vector<double> rho3(3);
		std::vector<double> rho(3);
		std::vector<double> rho0N(3);

		std::vector<double> KK(3);

		std::vector<double> rhotmp(3);
		std::vector<double> ll1(3);
		std::vector<double> ll2(3);
		std::vector<double> ll3(3);

		std::vector<double> u1(3);
		std::vector<double> u2(3);
		std::vector<double> u3(3);
		std::vector<double*> U(3);

		std::vector<double> vecpom(3);

		std::vector<double> RWGf(3);
		std::vector<double> RWGs(3);

		std::vector<double> p1nulvec(3);
		std::vector<double> p2nulvec(3);
		std::vector<double> p3nulvec(3);
		std::vector<double> pomVEC(3);
		std::vector<double> pomocvF(3);
		std::vector<double> pomocvS(3);
		std::vector<double> pomNvec(3);
		std::vector<double> pomVec1(3);
		std::vector<std::complex<double>> pomVec2(3);
		std::vector<std::complex<double>> pomVec3(3);
		std::vector<std::complex<double>> pomVec4(3);
		std::vector<double> pomVec5(3);
		std::vector<double> pomVec6(3);

		std::vector<double> Ivec1R(3);
		std::vector<double> Ivec1R3(3);
		std::vector<double> Rivec(3);
		std::vector<std::complex<double>> Ikon(3);
		std::vector<std::complex<double>> Ikon1R(3);
		std::vector<double> Ikon1R3(3);
		std::vector<std::complex<double>> KROSS1(3);
		std::vector<double> KROSS2(3);

		std::vector<double> POM(3);
		std::vector<double> POM1(3);
		std::vector<double> POM2(3);

		std::vector<double> Lpl(3);
		std::vector<double> Lmn(3);
		std::vector<double> Rpl(3);
		std::vector<double> Rmn(3);
		std::vector<double> Pnul(3);
		std::vector<double*> Pnulvec(3);
		std::vector<double> Rnul(3);
		std::vector<double> Ppl(3);
		std::vector<double> Pmn(3);
		std::vector<double> D(3);

		double const1, const2;
		std::complex<double> konst1, konst2;
		double LOGpl, LOGmn, ATANpl, ATANmn;
		double Isca1R, Isca1R3;
		double K;

		K = dot(&normalVectorS[0], &vertexS1[0]);
		mul(&npom1[0], &normalVectorS[0], K);
		sub(&rho1[0], &vertexS1[0], &npom1[0]);

		K = dot(&normalVectorS[0], &vertexS2[0]);
		mul(&npom2[0], &normalVectorS[0], K);
		sub(&rho2[0], &vertexS2[0], &npom2[0]);

		K = dot(&normalVectorS[0], &vertexS3[0]);
		mul(&npom3[0], &normalVectorS[0], K);
		sub(&rho3[0], &vertexS3[0], &npom3[0]);

		sub(&KK[0], &rho2[0], &rho1[0]);
		double norm1 = norm(&KK[0]);

		sub(&rhotmp[0], &rho2[0], &rho1[0]);
		div(&ll1[0], &rhotmp[0], norm1);
		sub(&rhotmp[0], &rho3[0], &rho2[0]);
		norm1 = norm(&rhotmp[0]);

		sub(&rhotmp[0], &rho3[0], &rho2[0]);
		div(&ll2[0], &rhotmp[0], norm1);
		sub(&rhotmp[0], &rho1[0], &rho3[0]);
		norm1 = norm(&rhotmp[0]);

		div(&ll3[0], &rhotmp[0], norm1);

		cross(&u1[0], &ll1[0], &normalVectorS[0]);
		cross(&u2[0], &ll2[0], &normalVectorS[0]);
		cross(&u3[0], &ll3[0], &normalVectorS[0]);

		U[0] = &u1[0];
		U[1] = &u2[0];
		U[2] = &u3[0];

		for (int numF = 0; numF < gp.N_points_; ++numF) {
			// ���㳡RWG��������������Ԫ�ϵĻ��ֵ�����
			fieldPoint[0] = gp.l1[numF] * vertexF1[0] + gp.l2[numF] * vertexF2[0] + gp.l3[numF] * vertexF3[0];
			fieldPoint[1] = gp.l1[numF] * vertexF1[1] + gp.l2[numF] * vertexF2[1] + gp.l3[numF] * vertexF3[1];
			fieldPoint[2] = gp.l1[numF] * vertexF1[2] + gp.l2[numF] * vertexF2[2] + gp.l3[numF] * vertexF3[2];
			// ��ȡ�û��ֵ�Ȩ��
			double wF = gp.weight[numF];

			K = dot(&normalVectorS[0], &fieldPoint[0]);
			mul(&vecpom[0], &normalVectorS[0], K);
			sub(&rho[0], &fieldPoint[0], &vecpom[0]);

			sub(&rhotmp[0], &rho2[0], &rho[0]);
			double l1pl = dot(&rhotmp[0], &ll1[0]);
			sub(&rhotmp[0], &rho1[0], &rho[0]);
			double l1mn = dot(&rhotmp[0], &ll1[0]);

			sub(&rhotmp[0], &rho3[0], &rho[0]);
			double l2pl = dot(&rhotmp[0], &ll2[0]);
			sub(&rhotmp[0], &rho2[0], &rho[0]);
			double l2mn = dot(&rhotmp[0], &ll2[0]);

			sub(&rhotmp[0], &rho1[0], &rho[0]);
			double l3pl = dot(&rhotmp[0], &ll3[0]);
			sub(&rhotmp[0], &rho3[0], &rho[0]);
			double l3mn = dot(&rhotmp[0], &ll3[0]);

			Lpl[0] = l1pl;
			Lpl[1] = l2pl;
			Lpl[2] = l3pl;
			Lmn[0] = l1mn;
			Lmn[1] = l2mn;
			Lmn[2] = l3mn;

			sub(&rhotmp[0], &rho2[0], &rho[0]);
			double p1nul = abs(dot(&rhotmp[0], &u1[0]));
			sub(&rhotmp[0], &rho3[0], &rho[0]);
			double p2nul = abs(dot(&rhotmp[0], &u2[0]));
			sub(&rhotmp[0], &rho1[0], &rho[0]);
			double p3nul = abs(dot(&rhotmp[0], &u3[0]));

			Pnul[0] = p1nul;
			Pnul[1] = p2nul;
			Pnul[2] = p3nul;

			sub(&rhotmp[0], &rho2[0], &rho[0]);
			double p1pl = norm(&rhotmp[0]);
			sub(&rhotmp[0], &rho1[0], &rho[0]);
			double p1mn = norm(&rhotmp[0]);

			sub(&rhotmp[0], &rho3[0], &rho[0]);
			double p2pl = norm(&rhotmp[0]);
			sub(&rhotmp[0], &rho2[0], &rho[0]);
			double p2mn = norm(&rhotmp[0]);

			sub(&rhotmp[0], &rho1[0], &rho[0]);
			double p3pl = norm(&rhotmp[0]);
			sub(&rhotmp[0], &rho3[0], &rho[0]);
			double p3mn = norm(&rhotmp[0]);

			Ppl[0] = p1pl;
			Ppl[1] = p2pl;
			Ppl[2] = p3pl;
			Pmn[0] = p1mn;
			Pmn[1] = p2mn;
			Pmn[2] = p3mn;

			if (p1nul < 50 * eps0) {
				p1nulvec[0] = 0.0;
				p1nulvec[1] = 0.0;
				p1nulvec[2] = 0.0;
			}
			else {
				sub(&POM1[0], &rho2[0], &rho[0]);
				mul(&POM2[0], &ll1[0], l1pl);
				sub(&POM[0], &POM1[0], &POM2[0]);
				div(&p1nulvec[0], &POM[0], p1nul);
			}

			if (p2nul < 50 * eps0) {
				p2nulvec[0] = 0.0;
				p2nulvec[1] = 0.0;
				p2nulvec[2] = 0.0;
			}
			else {
				sub(&POM1[0], &rho3[0], &rho[0]);
				mul(&POM2[0], &ll2[0], l2pl);
				sub(&POM[0], &POM1[0], &POM2[0]);
				div(&p2nulvec[0], &POM[0], p2nul);
			}

			if (p3nul < 50 * eps0) {
				p3nulvec[0] = 0.0;
				p3nulvec[1] = 0.0;
				p3nulvec[2] = 0.0;
			}
			else {
				sub(&POM1[0], &rho1[0], &rho[0]);
				mul(&POM2[0], &ll3[0], l3pl);
				sub(&POM[0], &POM1[0], &POM2[0]);
				div(&p3nulvec[0], &POM[0], p3nul);
			}

			Pnulvec[0] = &p1nulvec[0];
			Pnulvec[1] = &p2nulvec[0];
			Pnulvec[2] = &p3nulvec[0];

			sub(&rhotmp[0], &fieldPoint[0], vertexS1);
			double d1 = dot(&normalVectorS[0], rhotmp);
			sub(&rhotmp[0], &fieldPoint[0], vertexS2);
			double d2 = dot(&normalVectorS[0], rhotmp);
			sub(&rhotmp[0], &fieldPoint[0], vertexS3);
			double d3 = dot(&normalVectorS[0], rhotmp);

			D[0] = d1;
			D[1] = d2;
			D[2] = d3;

			double R1nul = sqrt(pow(p1nul, 2) + pow(d1, 2));
			double R2nul = sqrt(pow(p2nul, 2) + pow(d2, 2));
			double R3nul = sqrt(pow(p3nul, 2) + pow(d3, 2));

			Rnul[0] = R1nul;
			Rnul[1] = R2nul;
			Rnul[2] = R3nul;

			double R1pl = sqrt(pow(p1pl, 2) + pow(d1, 2));
			double R2pl = sqrt(pow(p2pl, 2) + pow(d2, 2));
			double R3pl = sqrt(pow(p3pl, 2) + pow(d3, 2));

			Rpl[0] = R1pl;
			Rpl[1] = R2pl;
			Rpl[2] = R3pl;

			double R1mn = sqrt(pow(p1mn, 2) + pow(d1, 2));
			double R2mn = sqrt(pow(p2mn, 2) + pow(d2, 2));
			double R3mn = sqrt(pow(p3mn, 2) + pow(d3, 2));

			Rmn[0] = R1mn;
			Rmn[1] = R2mn;
			Rmn[2] = R3mn;

			Ivec1R3[0] = Ivec1R3[1] = Ivec1R3[2] = Ivec1R[0] = Ivec1R[1] = Ivec1R[2] = 0.0;

			Isca1R = Isca1R3 = 0.0;

			// vec part 1/R singularity
			for (int iv = 0; iv < 3; ++iv) {
				if ((Rpl[iv] + Lpl[iv]) < 50 * eps0) {
					LOGpl = 0.0;
				}
				else {
					LOGpl = log(Rpl[iv] + Lpl[iv]);
				}

				if ((Rmn[iv] + Lmn[iv]) < 50 * eps0) {
					LOGmn = 0.0;
				}
				else {
					LOGmn = log(Rmn[iv] + Lmn[iv]);
				}

				const1 = 0.5 * (pow(Rnul[iv], 2) * (LOGpl - LOGmn) + Rpl[iv] * Lpl[iv] - Rmn[iv] * Lmn[iv]);

				mul(&pomVEC[0], U[iv], const1);
				add(&Ivec1R[0], &Ivec1R[0], &pomVEC[0]);
			}
			// scalar part 1/R singularity
			for (int is = 0; is < 3; ++is) {
				if ((Rpl[is] + Lpl[is]) < 50 * eps0) {
					LOGpl = 0.0;
				}
				else {
					LOGpl = log(Rpl[is] + Lpl[is]);
				}

				if ((Rmn[is] + Lmn[is]) < 50 * eps0) {
					LOGmn = 0.0;
				}
				else {
					LOGmn = log(Rmn[is] + Lmn[is]);
				}

				if ((pow(Rnul[is], 2) + abs(D[is]) * Rpl[is]) < 50 * eps0) {
					ATANpl = 0.0;
				}
				else {
					ATANpl = atan((Pnul[is] * Lpl[is]) / (pow(Rnul[is], 2) + abs(D[is]) * Rpl[is]));
				}

				if ((pow(Rnul[is], 2) + abs(D[is]) * Rmn[is]) < 50 * eps0) {
					ATANmn = 0.0;
				}
				else {
					ATANmn = atan((Pnul[is] * Lmn[is]) / (pow(Rnul[is], 2) + abs(D[is]) * Rmn[is]));
				}

				Isca1R = Isca1R + dot(Pnulvec[is], U[is]) * (Pnul[is] * (LOGpl - LOGmn) - abs(D[is]) * (ATANpl - ATANmn));
			}

			// vec part 1/R^3 singularity
			for (int iv = 0; iv < 3; ++iv) {
				if ((Rpl[iv] + Lpl[iv]) < 50.0 * eps0) {
					LOGpl = 0.0;
				}
				else {
					LOGpl = log(Rpl[iv] + Lpl[iv]);
				}

				if ((Rmn[iv] + Lmn[iv]) < 50.0 * eps0) {
					LOGmn = 0.0;
				}
				else {
					LOGmn = log(Rmn[iv] + Lmn[iv]);
				}

				const1 = LOGpl - LOGmn;

				mul(&pomVEC[0], U[iv], const1);
				sub(&Ivec1R3[0], &Ivec1R3[0], &pomVEC[0]);
			}

			// scalar part 1/R^3 singularity
			for (int is = 0; is < 3; ++is) {
				if (abs(D[is]) >= 10.0 * eps0) {
					ATANpl = atan((Pnul[is] * Lpl[is]) / (pow(Rnul[is], 2) + abs(D[is]) * Rpl[is]));
					ATANmn = atan((Pnul[is] * Lmn[is]) / (pow(Rnul[is], 2) + abs(D[is]) * Rmn[is]));
					Isca1R3 = Isca1R3 + dot(Pnulvec[is], U[is]) * ((1.0 / abs(D[is])) * (ATANpl - ATANmn));
				}
				else if ((abs(D[is]) < 10.0 * eps0) && (Pnul[is] >= 10.0 * eps0)) {
					Isca1R3 = Isca1R3 - dot(Pnulvec[is], U[is]) * (Lpl[is] / (Pnul[is] * Rpl[is]) - Lmn[is] / (Pnul[is] * Rmn[is]));
				}
				else {
					Isca1R3 = 0.0;
				}
			}

			// testing Galerkin
			sub(&pomocvF[0], &fieldPoint[0], nonPublicEdge_vertexF);

			mul(&RWGf[0], &pomocvF[0], recAf);

			sub(&Rivec[0], &fieldPoint[0], nonPublicEdge_vertexS);

			double K1 = dot(&normalVectorS[0], nonPublicEdge_vertexS);
			mul(&pomNvec[0], normalVectorS, K1);
			sub(&rho0N[0], nonPublicEdge_vertexS, &pomNvec[0]);

			konst1 = k_ * k_ / 2.0;
			konst2 = konst1 * Isca1R;

			sub(&pomVec1[0], &rho[0], &rho0N[0]);
			mul(&pomVec2[0], &Ivec1R[0], konst1);
			mul(&pomVec3[0], &pomVec1[0], konst2);
			add(&pomVec4[0], &pomVec2[0], &pomVec3[0]);

			mul(&Ikon1R[0], &pomVec4[0], recAs);

			mul(&pomVec5[0], &pomVec1[0], Isca1R3);
			add(&pomVec6[0], &Ivec1R3[0], &pomVec5[0]);
			mul(&Ikon1R3[0], &pomVec6[0], recAs);
			add(&Ikon[0], &Ikon1R3[0], &Ikon1R[0]);

			cross(&KROSS1[0], &Rivec[0], &Ikon[0]);
			cross(&KROSS2[0], &RWGf[0], normalVectorF);
			alok += wF * 2.0 * areaTriF * dot(&KROSS2[0], &KROSS1[0]);
		}
	}
}// end of namespace
```


## 8. K 算子的函数化提取

### 8.1 提取边界

MFIE 左端由两个不同数学对象组成：

```math
\frac12\mathbf J_s+\hat{\mathbf n}\times\mathcal K(\mathbf J_s).
```

因此函数化提取时既不能把 `1/2` 自作用项误认为 K 算子，也不能把 MFIE 外层的
`hat(n) ×` 混入 K 算子。这里严格按

```math
\mathcal K(\mathbf X)
=\iint_{S_0}\mathbf X(\mathbf r')\times
\nabla G_0(\mathbf r,\mathbf r')\,\mathrm dS'
```

提取函数。完整实现拆成三个函数：

- `computeK_Zij`：只计算 `K` 算子的 RWG 伽辽金矩阵元，不进行场面单位外法向量叉乘。
- `singularityK`：为 `computeK_Zij` 补回同面核中解析拆出的 `1/R^3` 与 `k^2/(2R)` 部分。
- `computeMFIEJump_Zij`：单独计算 MFIE 的 `1/2` 跳跃项。组装完整 MFIE 时，将它与 `computeK_Zij` 的结果相加。

| 原 MFIE 代码 | 函数化结果 |
|---|---|
| `computeMFIE_Zij` 中去掉外层 `hat(n) ×` 的 K 积分 | `KOperator::computeK_Zij` |
| `singularityMFIE` | `KOperator::singularityK` |
| 同面分支的单重 `numI` 积分 | `KOperator::computeMFIEJump_Zij` |
| `Zij_near_M` | `Zij_near_K` 或 `Zij_jump` |

`computeK_Zij` 直接以场 RWG 测试函数点乘
`R × f_n`，不再计算 `normalVectorF × (R × f_n)`。
`singularityK` 同样直接测试解析补偿所得的 K 向量，不再执行场面外法向量叉乘。
其中保留的 `normalVectorS` 仅用于建立源三角形的局部几何坐标和求解析奇异积分，
并不参与 `hat(n) × K` 运算。

### 8.2 完整提取头文件

```cpp
// KOperator.h
#pragma once
#ifndef K_OPERATOR_H
#define K_OPERATOR_H

#include "MyStruct.h"
#include "GaussPoints.h"

namespace KOperator {
	// 仅计算 K 算子矩阵元：不包含 MFIE 的 1/2 自作用项，
	// 也不包含场面单位外法向量与 K 的叉乘。
	void computeK_Zij(const RWGBase* rwgField, const RWGBase* rwgSource,
		const std::complex<double> k_, const gaussPoints& gp,
		std::complex<double>& Zij_near_K);

	// 单独计算 MFIE 的 1/2 电流跳跃项。
	void computeMFIEJump_Zij(const RWGBase* rwgField, const RWGBase* rwgSource,
		const gaussPoints& gp, std::complex<double>& Zij_jump);

	// 计算 K 算子同面 1/R^3 与 1/R 奇异核的解析补偿。
	void singularityK(double& areaTriF, double& areaTriS,
		double* vertexF1, double* vertexF2, double* vertexF3,
		double* nonPublicEdge_vertexF, const double recAf,
		double* vertexS1, double* vertexS2, double* vertexS3,
		double* nonPublicEdge_vertexS, double* normalVectorS, const double recAs,
		const gaussPoints& gp, const std::complex<double> k_,
		std::complex<double>& alok);
}
#endif // !K_OPERATOR_H
```

### 8.3 完整提取实现文件

```cpp
// KOperator.cpp
#include "KOperator.h"
#include "OverloadAlgo.h"

namespace KOperator {
	void computeK_Zij(const RWGBase* rwgField, const RWGBase* rwgSource, const std::complex<double> k_,
		const gaussPoints& gp, std::complex<double>& Zij_near_K)
	{
		int rwgFieldTriId, rwgSourceTriId;
		double areaTriF, areaTriS;
		double rec2Af, rec2As;
		double* nonPublicEdge_vertexF, * nonPublicEdge_vertexS;
		double* vertexF1, * vertexF2, * vertexF3;
		double* vertexS1, * vertexS2, * vertexS3;
		double* normalVectorS;
		double lF, lS;

		std::vector<double> fieldPoint(3), pomocvF(3), RWGf(3);
		std::vector<double> sourcePoint(3), pomocvS(3), RWGs(3);

		std::vector<double> CROSS1(3);

		std::vector<double> Rvec(3);
		std::vector<double> Rivec(3);

		for (int fi = 0; fi < 2; ++fi) {
			lF = rwgField->edgeLength;
			if (fi == 0) {
				rwgFieldTriId = rwgField->negativeFace->faceId;
				// ��RWG�����������������
				areaTriF = rwgField->negativeFace->triangularArea;
				rec2Af = 1.0 / (2.0 * areaTriF);
				// ע�⣺����ķǹ����߶����Ǹ������εĶ���
				nonPublicEdge_vertexF = &(rwgField->freeVertexNegative->x);
				// ��ó�RWG������������������
				vertexF1 = &(rwgField->negativeFace->vertex1->x);
				vertexF2 = &(rwgField->negativeFace->vertex2->x);
				vertexF3 = &(rwgField->negativeFace->vertex3->x);
			}
			else {
				rwgFieldTriId = rwgField->positiveFace->faceId;
				// ��RWG�����������������
				areaTriF = rwgField->positiveFace->triangularArea;
				rec2Af = 1.0 / (2.0 * areaTriF);
				// ��ȡ��RWG�������ķǹ����߶���
				nonPublicEdge_vertexF = &(rwgField->freeVertexPositive->x);
				// ��ó�RWG������������������
				vertexF1 = &(rwgField->positiveFace->vertex1->x);
				vertexF2 = &(rwgField->positiveFace->vertex2->x);
				vertexF3 = &(rwgField->positiveFace->vertex3->x);
			}

			for (int sj = 0; sj < 2; ++sj) {
				lS = rwgSource->edgeLength;
				if (sj == 0) {
					rwgSourceTriId = rwgSource->negativeFace->faceId;
					// ԴRWG�����������������
					areaTriS = rwgSource->negativeFace->triangularArea;
					rec2As = 1.0 / (2.0 * areaTriS);
					// ��ȡ��RWG�������ķǹ����߶���
					// ע�⣺����ķǹ����߶����Ǹ������εĶ���
					nonPublicEdge_vertexS = &(rwgSource->freeVertexNegative->x);
					// ��ó�RWG������������������
					vertexS1 = &(rwgSource->negativeFace->vertex1->x);
					vertexS2 = &(rwgSource->negativeFace->vertex2->x);
					vertexS3 = &(rwgSource->negativeFace->vertex3->x);
					// ���ԴRWG�������������εķ�����
					normalVectorS = &(rwgSource->negativeFace->externalNormalVector.x);
				}
				else {
					rwgSourceTriId = rwgSource->positiveFace->faceId;
					//  ԴRWG�����������������
					areaTriS = rwgSource->positiveFace->triangularArea;
					rec2As = 1.0 / (2.0 * areaTriS);
					// ��ȡԴRWG�������ķǹ����߶���
					nonPublicEdge_vertexS = &(rwgSource->freeVertexPositive->x);
					// ���ԴRWG������������������
					vertexS1 = &(rwgSource->positiveFace->vertex1->x);
					vertexS2 = &(rwgSource->positiveFace->vertex2->x);
					vertexS3 = &(rwgSource->positiveFace->vertex3->x);
					// ���ԴRWG�������������εķ�����
					normalVectorS = &(rwgSource->positiveFace->externalNormalVector.x);
				}

				std::complex<double> alok(0.0, 0.0);
				std::complex<double> alok1(0.0, 0.0);

				if (rwgFieldTriId != rwgSourceTriId) {
					// ���ù�ʽ: $\sum_{f=1}^{n_P}\sum_{s=1}^{n_P}w_fw_s\mathbf{f}_m(\mathbf{r}_f)\cdot\left[\hat{n}_f\times((\mathbf{r}_f-\mathbf{r}_s)/R^3(1+ikR)e^{-ikR}\times\mathbf{f}_n(\mathbf{r}_s))\right]$
					// ѭ����������������Ԫ�ϵĻ��ֵ�
					for (int numI = 0; numI < gp.N_points_; ++numI) {
						// ���㳡RWG��������������Ԫ�ϵĻ��ֵ�����
						fieldPoint[0] = gp.l1[numI] * vertexF1[0] + gp.l2[numI] * vertexF2[0] + gp.l3[numI] * vertexF3[0];
						fieldPoint[1] = gp.l1[numI] * vertexF1[1] + gp.l2[numI] * vertexF2[1] + gp.l3[numI] * vertexF3[1];
						fieldPoint[2] = gp.l1[numI] * vertexF1[2] + gp.l2[numI] * vertexF2[2] + gp.l3[numI] * vertexF3[2];
						// ��ȡ�ó����ֵ�Ȩ��
						double wF = gp.weight[numI];

						// ѭ������Դ��������Ԫ�ϵĻ��ֵ�
						for (int numJ = 0; numJ < gp.N_points_; ++numJ) {
							// ����ԴRWG��������������Ԫ�ϵĻ��ֵ�����
							sourcePoint[0] = gp.l1[numJ] * vertexS1[0] + gp.l2[numJ] * vertexS2[0] + gp.l3[numJ] * vertexS3[0];
							sourcePoint[1] = gp.l1[numJ] * vertexS1[1] + gp.l2[numJ] * vertexS2[1] + gp.l3[numJ] * vertexS3[1];
							sourcePoint[2] = gp.l1[numJ] * vertexS1[2] + gp.l2[numJ] * vertexS2[2] + gp.l3[numJ] * vertexS3[2];
							// ��ȡ��Դ���ֵ�Ȩ��
							double wS = gp.weight[numJ];
							// ����ѡ����RWG��������������Ԫ�ϵĻ��ֵ㵽ԴRWG��������������Ԫ�ϵĻ��ֵ�ľ���
							// ��˹���ֹ���߻�ѡ��ͬһ��
							double R = norm(&fieldPoint[0], &sourcePoint[0]);

							sub(&Rvec[0], &fieldPoint[0], &sourcePoint[0]);

							std::complex<double> gradG = (exp(-J * k_ * R) / (R * R * R)) * (1.0 + J * k_ * R);

							sub(&pomocvF[0], &fieldPoint[0], nonPublicEdge_vertexF);
							sub(&pomocvS[0], &sourcePoint[0], nonPublicEdge_vertexS);

							cross(&CROSS1[0], &Rvec[0], &pomocvS[0]);
							std::complex<double> contrib = wF * wS * (gradG * dot(&pomocvF[0], &CROSS1[0]));
							alok1 += contrib;
						}
					}
				}
				else {
					// ���ù�ʽ: $\sum_{f}\sum_{f}4w_{s}w_{s}A_{f}A_{s}\left[\nabla G(\mathbf{r}_{f},\mathbf{r}_{s})\cdot\left((\mathbf{r}_{f}-\mathbf{r}_{s})\times\mathbf{f}_{n}(\mathbf{r}_{s})\cdot(\mathbf{f}_{m}(\mathbf{r}_{f})\times\hat{n}_{f})\right)\right]$
					for (int numI = 0; numI < gp.N_points_; ++numI) {
						// ���㳡RWG��������������Ԫ�ϵĻ��ֵ�����
						fieldPoint[0] = gp.l1[numI] * vertexF1[0] + gp.l2[numI] * vertexF2[0] + gp.l3[numI] * vertexF3[0];
						fieldPoint[1] = gp.l1[numI] * vertexF1[1] + gp.l2[numI] * vertexF2[1] + gp.l3[numI] * vertexF3[1];
						fieldPoint[2] = gp.l1[numI] * vertexF1[2] + gp.l2[numI] * vertexF2[2] + gp.l3[numI] * vertexF3[2];
						// ��ȡ�ó����ֵ�Ȩ��
						double wF = gp.weight[numI];

						// ѭ������Դ��������Ԫ�ϵĻ��ֵ�
						for (int numJ = 0; numJ < gp.N_points_; ++numJ) {
							// ����ԴRWG��������������Ԫ�ϵĻ��ֵ�����
							sourcePoint[0] = gp.l1[numJ] * vertexS1[0] + gp.l2[numJ] * vertexS2[0] + gp.l3[numJ] * vertexS3[0];
							sourcePoint[1] = gp.l1[numJ] * vertexS1[1] + gp.l2[numJ] * vertexS2[1] + gp.l3[numJ] * vertexS3[1];
							sourcePoint[2] = gp.l1[numJ] * vertexS1[2] + gp.l2[numJ] * vertexS2[2] + gp.l3[numJ] * vertexS3[2];
							// ��ȡ��Դ���ֵ�Ȩ��
							double wS = gp.weight[numJ];

							// ����ѡ����RWG��������������Ԫ�ϵĻ��ֵ㵽ԴRWG��������������Ԫ�ϵĻ��ֵ�ľ���
							double R = norm(&fieldPoint[0], &sourcePoint[0]);

							sub(&Rvec[0], &fieldPoint[0], &sourcePoint[0]);

							std::complex<double> G;

							if (R < eps0) {
								G = std::complex<double>(0.0, 0.0);
							}
							else {
								G = (exp(-J * k_ * R) / (R * R * R)) * (1.0 + J * k_ * R) - 1.0 / (R * R * R) - 0.5 * (k_ * k_) * (1.0 / R);
							}

							sub(&pomocvF[0], &fieldPoint[0], nonPublicEdge_vertexF);
							sub(&pomocvS[0], &sourcePoint[0], nonPublicEdge_vertexS);

							sub(&Rivec[0], &fieldPoint[0], nonPublicEdge_vertexS);

							cross(&CROSS1[0], &Rivec[0], &pomocvS[0]);
							std::complex<double> contrib = wF * wS * (G * dot(&pomocvF[0], &CROSS1[0]));
							alok1 += contrib;
							//temp2 += contrib;
						}
					}

					singularityK(areaTriF, areaTriS,
						vertexF1, vertexF2, vertexF3,
						nonPublicEdge_vertexF, rec2Af,
						vertexS1, vertexS2, vertexS3,
						nonPublicEdge_vertexS, normalVectorS, rec2As,
						gp, k_, alok);

					alok1 += alok;
				}

				// ��Ԫ�������������������Σ���Ϊ1������Ϊ-1
				const double signSymbolF = (fi == 0) ? -1.0 : 1.0;
				const double signSymbolS = (sj == 0) ? -1.0 : 1.0;
				Zij_near_K += (lF * lS) / (4.0 * Pi) * signSymbolF * signSymbolS * alok1;
			}
		}
	}

	void computeMFIEJump_Zij(const RWGBase* rwgField, const RWGBase* rwgSource,
		const gaussPoints& gp, std::complex<double>& Zij_jump)
	{
		std::vector<double> fieldPoint(3);
		std::vector<double> rhoField(3), rhoSource(3);
		std::vector<double> rwgFieldValue(3), rwgSourceValue(3);

		for (int fi = 0; fi < 2; ++fi) {
			auto* faceF = (fi == 0) ? rwgField->negativeFace : rwgField->positiveFace;
			double* freeVertexF = (fi == 0)
				? &(rwgField->freeVertexNegative->x)
				: &(rwgField->freeVertexPositive->x);
			double* vertexF1 = &(faceF->vertex1->x);
			double* vertexF2 = &(faceF->vertex2->x);
			double* vertexF3 = &(faceF->vertex3->x);
			const double areaF = faceF->triangularArea;
			const double rec2Af = 1.0 / (2.0 * areaF);

			for (int sj = 0; sj < 2; ++sj) {
				auto* faceS = (sj == 0) ? rwgSource->negativeFace : rwgSource->positiveFace;
				if (faceF->faceId != faceS->faceId) {
					continue;
				}

				double* freeVertexS = (sj == 0)
					? &(rwgSource->freeVertexNegative->x)
					: &(rwgSource->freeVertexPositive->x);
				const double areaS = faceS->triangularArea;
				const double rec2As = 1.0 / (2.0 * areaS);
				std::complex<double> localJump(0.0, 0.0);

				for (int numI = 0; numI < gp.N_points_; ++numI) {
					fieldPoint[0] = gp.l1[numI] * vertexF1[0]
						+ gp.l2[numI] * vertexF2[0]
						+ gp.l3[numI] * vertexF3[0];
					fieldPoint[1] = gp.l1[numI] * vertexF1[1]
						+ gp.l2[numI] * vertexF2[1]
						+ gp.l3[numI] * vertexF3[1];
					fieldPoint[2] = gp.l1[numI] * vertexF1[2]
						+ gp.l2[numI] * vertexF2[2]
						+ gp.l3[numI] * vertexF3[2];

					sub(&rhoField[0], &fieldPoint[0], freeVertexF);
					mul(&rwgFieldValue[0], &rhoField[0], rec2Af);
					sub(&rhoSource[0], &fieldPoint[0], freeVertexS);
					mul(&rwgSourceValue[0], &rhoSource[0], rec2As);

					localJump += 4.0 * Pi * gp.weight[numI] * areaF
						* dot(&rwgFieldValue[0], &rwgSourceValue[0]);
				}

				const double signF = (fi == 0) ? -1.0 : 1.0;
				const double signS = (sj == 0) ? -1.0 : 1.0;
				Zij_jump += rwgField->edgeLength * rwgSource->edgeLength
					/ (4.0 * Pi) * signF * signS * localJump;
			}
		}
	}

	void singularityK(double& areaTriF, double& areaTriS,
		double* vertexF1, double* vertexF2, double* vertexF3,
		double* nonPublicEdge_vertexF, const double recAf,
		double* vertexS1, double* vertexS2, double* vertexS3,
		double* nonPublicEdge_vertexS, double* normalVectorS, const double recAs,
		const gaussPoints& gp, const std::complex<double> k_, std::complex<double>& alok)
	{
		std::vector<double> fieldPoint(3);

		std::vector<double> npom1(3);
		std::vector<double> npom2(3);
		std::vector<double> npom3(3);
		std::vector<double> rho1(3);
		std::vector<double> rho2(3);
		std::vector<double> rho3(3);
		std::vector<double> rho(3);
		std::vector<double> rho0N(3);

		std::vector<double> KK(3);

		std::vector<double> rhotmp(3);
		std::vector<double> ll1(3);
		std::vector<double> ll2(3);
		std::vector<double> ll3(3);

		std::vector<double> u1(3);
		std::vector<double> u2(3);
		std::vector<double> u3(3);
		std::vector<double*> U(3);

		std::vector<double> vecpom(3);

		std::vector<double> RWGf(3);
		std::vector<double> RWGs(3);

		std::vector<double> p1nulvec(3);
		std::vector<double> p2nulvec(3);
		std::vector<double> p3nulvec(3);
		std::vector<double> pomVEC(3);
		std::vector<double> pomocvF(3);
		std::vector<double> pomocvS(3);
		std::vector<double> pomNvec(3);
		std::vector<double> pomVec1(3);
		std::vector<std::complex<double>> pomVec2(3);
		std::vector<std::complex<double>> pomVec3(3);
		std::vector<std::complex<double>> pomVec4(3);
		std::vector<double> pomVec5(3);
		std::vector<double> pomVec6(3);

		std::vector<double> Ivec1R(3);
		std::vector<double> Ivec1R3(3);
		std::vector<double> Rivec(3);
		std::vector<std::complex<double>> Ikon(3);
		std::vector<std::complex<double>> Ikon1R(3);
		std::vector<double> Ikon1R3(3);
		std::vector<std::complex<double>> KROSS1(3);

		std::vector<double> POM(3);
		std::vector<double> POM1(3);
		std::vector<double> POM2(3);

		std::vector<double> Lpl(3);
		std::vector<double> Lmn(3);
		std::vector<double> Rpl(3);
		std::vector<double> Rmn(3);
		std::vector<double> Pnul(3);
		std::vector<double*> Pnulvec(3);
		std::vector<double> Rnul(3);
		std::vector<double> Ppl(3);
		std::vector<double> Pmn(3);
		std::vector<double> D(3);

		double const1, const2;
		std::complex<double> konst1, konst2;
		double LOGpl, LOGmn, ATANpl, ATANmn;
		double Isca1R, Isca1R3;
		double K;

		K = dot(&normalVectorS[0], &vertexS1[0]);
		mul(&npom1[0], &normalVectorS[0], K);
		sub(&rho1[0], &vertexS1[0], &npom1[0]);

		K = dot(&normalVectorS[0], &vertexS2[0]);
		mul(&npom2[0], &normalVectorS[0], K);
		sub(&rho2[0], &vertexS2[0], &npom2[0]);

		K = dot(&normalVectorS[0], &vertexS3[0]);
		mul(&npom3[0], &normalVectorS[0], K);
		sub(&rho3[0], &vertexS3[0], &npom3[0]);

		sub(&KK[0], &rho2[0], &rho1[0]);
		double norm1 = norm(&KK[0]);

		sub(&rhotmp[0], &rho2[0], &rho1[0]);
		div(&ll1[0], &rhotmp[0], norm1);
		sub(&rhotmp[0], &rho3[0], &rho2[0]);
		norm1 = norm(&rhotmp[0]);

		sub(&rhotmp[0], &rho3[0], &rho2[0]);
		div(&ll2[0], &rhotmp[0], norm1);
		sub(&rhotmp[0], &rho1[0], &rho3[0]);
		norm1 = norm(&rhotmp[0]);

		div(&ll3[0], &rhotmp[0], norm1);

		cross(&u1[0], &ll1[0], &normalVectorS[0]);
		cross(&u2[0], &ll2[0], &normalVectorS[0]);
		cross(&u3[0], &ll3[0], &normalVectorS[0]);

		U[0] = &u1[0];
		U[1] = &u2[0];
		U[2] = &u3[0];

		for (int numF = 0; numF < gp.N_points_; ++numF) {
			// ���㳡RWG��������������Ԫ�ϵĻ��ֵ�����
			fieldPoint[0] = gp.l1[numF] * vertexF1[0] + gp.l2[numF] * vertexF2[0] + gp.l3[numF] * vertexF3[0];
			fieldPoint[1] = gp.l1[numF] * vertexF1[1] + gp.l2[numF] * vertexF2[1] + gp.l3[numF] * vertexF3[1];
			fieldPoint[2] = gp.l1[numF] * vertexF1[2] + gp.l2[numF] * vertexF2[2] + gp.l3[numF] * vertexF3[2];
			// ��ȡ�û��ֵ�Ȩ��
			double wF = gp.weight[numF];

			K = dot(&normalVectorS[0], &fieldPoint[0]);
			mul(&vecpom[0], &normalVectorS[0], K);
			sub(&rho[0], &fieldPoint[0], &vecpom[0]);

			sub(&rhotmp[0], &rho2[0], &rho[0]);
			double l1pl = dot(&rhotmp[0], &ll1[0]);
			sub(&rhotmp[0], &rho1[0], &rho[0]);
			double l1mn = dot(&rhotmp[0], &ll1[0]);

			sub(&rhotmp[0], &rho3[0], &rho[0]);
			double l2pl = dot(&rhotmp[0], &ll2[0]);
			sub(&rhotmp[0], &rho2[0], &rho[0]);
			double l2mn = dot(&rhotmp[0], &ll2[0]);

			sub(&rhotmp[0], &rho1[0], &rho[0]);
			double l3pl = dot(&rhotmp[0], &ll3[0]);
			sub(&rhotmp[0], &rho3[0], &rho[0]);
			double l3mn = dot(&rhotmp[0], &ll3[0]);

			Lpl[0] = l1pl;
			Lpl[1] = l2pl;
			Lpl[2] = l3pl;
			Lmn[0] = l1mn;
			Lmn[1] = l2mn;
			Lmn[2] = l3mn;

			sub(&rhotmp[0], &rho2[0], &rho[0]);
			double p1nul = abs(dot(&rhotmp[0], &u1[0]));
			sub(&rhotmp[0], &rho3[0], &rho[0]);
			double p2nul = abs(dot(&rhotmp[0], &u2[0]));
			sub(&rhotmp[0], &rho1[0], &rho[0]);
			double p3nul = abs(dot(&rhotmp[0], &u3[0]));

			Pnul[0] = p1nul;
			Pnul[1] = p2nul;
			Pnul[2] = p3nul;

			sub(&rhotmp[0], &rho2[0], &rho[0]);
			double p1pl = norm(&rhotmp[0]);
			sub(&rhotmp[0], &rho1[0], &rho[0]);
			double p1mn = norm(&rhotmp[0]);

			sub(&rhotmp[0], &rho3[0], &rho[0]);
			double p2pl = norm(&rhotmp[0]);
			sub(&rhotmp[0], &rho2[0], &rho[0]);
			double p2mn = norm(&rhotmp[0]);

			sub(&rhotmp[0], &rho1[0], &rho[0]);
			double p3pl = norm(&rhotmp[0]);
			sub(&rhotmp[0], &rho3[0], &rho[0]);
			double p3mn = norm(&rhotmp[0]);

			Ppl[0] = p1pl;
			Ppl[1] = p2pl;
			Ppl[2] = p3pl;
			Pmn[0] = p1mn;
			Pmn[1] = p2mn;
			Pmn[2] = p3mn;

			if (p1nul < 50 * eps0) {
				p1nulvec[0] = 0.0;
				p1nulvec[1] = 0.0;
				p1nulvec[2] = 0.0;
			}
			else {
				sub(&POM1[0], &rho2[0], &rho[0]);
				mul(&POM2[0], &ll1[0], l1pl);
				sub(&POM[0], &POM1[0], &POM2[0]);
				div(&p1nulvec[0], &POM[0], p1nul);
			}

			if (p2nul < 50 * eps0) {
				p2nulvec[0] = 0.0;
				p2nulvec[1] = 0.0;
				p2nulvec[2] = 0.0;
			}
			else {
				sub(&POM1[0], &rho3[0], &rho[0]);
				mul(&POM2[0], &ll2[0], l2pl);
				sub(&POM[0], &POM1[0], &POM2[0]);
				div(&p2nulvec[0], &POM[0], p2nul);
			}

			if (p3nul < 50 * eps0) {
				p3nulvec[0] = 0.0;
				p3nulvec[1] = 0.0;
				p3nulvec[2] = 0.0;
			}
			else {
				sub(&POM1[0], &rho1[0], &rho[0]);
				mul(&POM2[0], &ll3[0], l3pl);
				sub(&POM[0], &POM1[0], &POM2[0]);
				div(&p3nulvec[0], &POM[0], p3nul);
			}

			Pnulvec[0] = &p1nulvec[0];
			Pnulvec[1] = &p2nulvec[0];
			Pnulvec[2] = &p3nulvec[0];

			sub(&rhotmp[0], &fieldPoint[0], vertexS1);
			double d1 = dot(&normalVectorS[0], rhotmp);
			sub(&rhotmp[0], &fieldPoint[0], vertexS2);
			double d2 = dot(&normalVectorS[0], rhotmp);
			sub(&rhotmp[0], &fieldPoint[0], vertexS3);
			double d3 = dot(&normalVectorS[0], rhotmp);

			D[0] = d1;
			D[1] = d2;
			D[2] = d3;

			double R1nul = sqrt(pow(p1nul, 2) + pow(d1, 2));
			double R2nul = sqrt(pow(p2nul, 2) + pow(d2, 2));
			double R3nul = sqrt(pow(p3nul, 2) + pow(d3, 2));

			Rnul[0] = R1nul;
			Rnul[1] = R2nul;
			Rnul[2] = R3nul;

			double R1pl = sqrt(pow(p1pl, 2) + pow(d1, 2));
			double R2pl = sqrt(pow(p2pl, 2) + pow(d2, 2));
			double R3pl = sqrt(pow(p3pl, 2) + pow(d3, 2));

			Rpl[0] = R1pl;
			Rpl[1] = R2pl;
			Rpl[2] = R3pl;

			double R1mn = sqrt(pow(p1mn, 2) + pow(d1, 2));
			double R2mn = sqrt(pow(p2mn, 2) + pow(d2, 2));
			double R3mn = sqrt(pow(p3mn, 2) + pow(d3, 2));

			Rmn[0] = R1mn;
			Rmn[1] = R2mn;
			Rmn[2] = R3mn;

			Ivec1R3[0] = Ivec1R3[1] = Ivec1R3[2] = Ivec1R[0] = Ivec1R[1] = Ivec1R[2] = 0.0;

			Isca1R = Isca1R3 = 0.0;

			// vec part 1/R singularity
			for (int iv = 0; iv < 3; ++iv) {
				if ((Rpl[iv] + Lpl[iv]) < 50 * eps0) {
					LOGpl = 0.0;
				}
				else {
					LOGpl = log(Rpl[iv] + Lpl[iv]);
				}

				if ((Rmn[iv] + Lmn[iv]) < 50 * eps0) {
					LOGmn = 0.0;
				}
				else {
					LOGmn = log(Rmn[iv] + Lmn[iv]);
				}

				const1 = 0.5 * (pow(Rnul[iv], 2) * (LOGpl - LOGmn) + Rpl[iv] * Lpl[iv] - Rmn[iv] * Lmn[iv]);

				mul(&pomVEC[0], U[iv], const1);
				add(&Ivec1R[0], &Ivec1R[0], &pomVEC[0]);
			}
			// scalar part 1/R singularity
			for (int is = 0; is < 3; ++is) {
				if ((Rpl[is] + Lpl[is]) < 50 * eps0) {
					LOGpl = 0.0;
				}
				else {
					LOGpl = log(Rpl[is] + Lpl[is]);
				}

				if ((Rmn[is] + Lmn[is]) < 50 * eps0) {
					LOGmn = 0.0;
				}
				else {
					LOGmn = log(Rmn[is] + Lmn[is]);
				}

				if ((pow(Rnul[is], 2) + abs(D[is]) * Rpl[is]) < 50 * eps0) {
					ATANpl = 0.0;
				}
				else {
					ATANpl = atan((Pnul[is] * Lpl[is]) / (pow(Rnul[is], 2) + abs(D[is]) * Rpl[is]));
				}

				if ((pow(Rnul[is], 2) + abs(D[is]) * Rmn[is]) < 50 * eps0) {
					ATANmn = 0.0;
				}
				else {
					ATANmn = atan((Pnul[is] * Lmn[is]) / (pow(Rnul[is], 2) + abs(D[is]) * Rmn[is]));
				}

				Isca1R = Isca1R + dot(Pnulvec[is], U[is]) * (Pnul[is] * (LOGpl - LOGmn) - abs(D[is]) * (ATANpl - ATANmn));
			}

			// vec part 1/R^3 singularity
			for (int iv = 0; iv < 3; ++iv) {
				if ((Rpl[iv] + Lpl[iv]) < 50.0 * eps0) {
					LOGpl = 0.0;
				}
				else {
					LOGpl = log(Rpl[iv] + Lpl[iv]);
				}

				if ((Rmn[iv] + Lmn[iv]) < 50.0 * eps0) {
					LOGmn = 0.0;
				}
				else {
					LOGmn = log(Rmn[iv] + Lmn[iv]);
				}

				const1 = LOGpl - LOGmn;

				mul(&pomVEC[0], U[iv], const1);
				sub(&Ivec1R3[0], &Ivec1R3[0], &pomVEC[0]);
			}

			// scalar part 1/R^3 singularity
			for (int is = 0; is < 3; ++is) {
				if (abs(D[is]) >= 10.0 * eps0) {
					ATANpl = atan((Pnul[is] * Lpl[is]) / (pow(Rnul[is], 2) + abs(D[is]) * Rpl[is]));
					ATANmn = atan((Pnul[is] * Lmn[is]) / (pow(Rnul[is], 2) + abs(D[is]) * Rmn[is]));
					Isca1R3 = Isca1R3 + dot(Pnulvec[is], U[is]) * ((1.0 / abs(D[is])) * (ATANpl - ATANmn));
				}
				else if ((abs(D[is]) < 10.0 * eps0) && (Pnul[is] >= 10.0 * eps0)) {
					Isca1R3 = Isca1R3 - dot(Pnulvec[is], U[is]) * (Lpl[is] / (Pnul[is] * Rpl[is]) - Lmn[is] / (Pnul[is] * Rmn[is]));
				}
				else {
					Isca1R3 = 0.0;
				}
			}

			// testing Galerkin
			sub(&pomocvF[0], &fieldPoint[0], nonPublicEdge_vertexF);

			mul(&RWGf[0], &pomocvF[0], recAf);

			sub(&Rivec[0], &fieldPoint[0], nonPublicEdge_vertexS);

			double K1 = dot(&normalVectorS[0], nonPublicEdge_vertexS);
			mul(&pomNvec[0], normalVectorS, K1);
			sub(&rho0N[0], nonPublicEdge_vertexS, &pomNvec[0]);

			konst1 = k_ * k_ / 2.0;
			konst2 = konst1 * Isca1R;

			sub(&pomVec1[0], &rho[0], &rho0N[0]);
			mul(&pomVec2[0], &Ivec1R[0], konst1);
			mul(&pomVec3[0], &pomVec1[0], konst2);
			add(&pomVec4[0], &pomVec2[0], &pomVec3[0]);

			mul(&Ikon1R[0], &pomVec4[0], recAs);

			mul(&pomVec5[0], &pomVec1[0], Isca1R3);
			add(&pomVec6[0], &Ivec1R3[0], &pomVec5[0]);
			mul(&Ikon1R3[0], &pomVec6[0], recAs);
			add(&Ikon[0], &Ikon1R3[0], &Ikon1R[0]);

			cross(&KROSS1[0], &Rivec[0], &Ikon[0]);
			alok += wF * 2.0 * areaTriF * dot(&RWGf[0], &KROSS1[0]);
		}
	}
}// end of namespace
```
