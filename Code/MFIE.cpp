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
				// 场RWG基函数负三角形面积
				areaTriF = rwgField->negativeFace->triangularArea;
				rec2Af = 1.0 / (2.0 * areaTriF);
				// 注意：这里的非公共边顶点是负三角形的顶点
				nonPublicEdge_vertexF = &(rwgField->freeVertexNegative->x);
				// 获得场RWG基函数的三顶点坐标
				vertexF1 = &(rwgField->negativeFace->vertex1->x);
				vertexF2 = &(rwgField->negativeFace->vertex2->x);
				vertexF3 = &(rwgField->negativeFace->vertex3->x);
				normalVectorF = &(rwgField->negativeFace->externalNormalVector.x);
			}
			else {
				rwgFieldTriId = rwgField->positiveFace->faceId;
				// 场RWG基函数正三角形面积
				areaTriF = rwgField->positiveFace->triangularArea;
				rec2Af = 1.0 / (2.0 * areaTriF);
				// 获取场RWG基函数的非公共边顶点
				nonPublicEdge_vertexF = &(rwgField->freeVertexPositive->x);
				// 获得场RWG基函数的三顶点坐标
				vertexF1 = &(rwgField->positiveFace->vertex1->x);
				vertexF2 = &(rwgField->positiveFace->vertex2->x);
				vertexF3 = &(rwgField->positiveFace->vertex3->x);
				normalVectorF = &(rwgField->positiveFace->externalNormalVector.x);
			}

			for (int sj = 0; sj < 2; ++sj) {
				lS = rwgSource->edgeLength;
				if (sj == 0) {
					rwgSourceTriId = rwgSource->negativeFace->faceId;
					// 源RWG基函数负三角形面积
					areaTriS = rwgSource->negativeFace->triangularArea;
					rec2As = 1.0 / (2.0 * areaTriS);
					// 获取场RWG基函数的非公共边顶点
					// 注意：这里的非公共边顶点是负三角形的顶点
					nonPublicEdge_vertexS = &(rwgSource->freeVertexNegative->x);
					// 获得场RWG基函数的三顶点坐标
					vertexS1 = &(rwgSource->negativeFace->vertex1->x);
					vertexS2 = &(rwgSource->negativeFace->vertex2->x);
					vertexS3 = &(rwgSource->negativeFace->vertex3->x);
					// 获得源RWG基函数负三角形的法向量
					normalVectorS = &(rwgSource->negativeFace->externalNormalVector.x);
				}
				else {
					rwgSourceTriId = rwgSource->positiveFace->faceId;
					//  源RWG基函数正三角形面积
					areaTriS = rwgSource->positiveFace->triangularArea;
					rec2As = 1.0 / (2.0 * areaTriS);
					// 获取源RWG基函数的非公共边顶点
					nonPublicEdge_vertexS = &(rwgSource->freeVertexPositive->x);
					// 获得源RWG基函数的三顶点坐标
					vertexS1 = &(rwgSource->positiveFace->vertex1->x);
					vertexS2 = &(rwgSource->positiveFace->vertex2->x);
					vertexS3 = &(rwgSource->positiveFace->vertex3->x);
					// 获得源RWG基函数正三角形的法向量
					normalVectorS = &(rwgSource->positiveFace->externalNormalVector.x);
				}

				std::complex<double> alok(0.0, 0.0);
				std::complex<double> alok1(0.0, 0.0);

				if (rwgFieldTriId != rwgSourceTriId) {
					// 所用公式: $\sum_{f=1}^{n_P}\sum_{s=1}^{n_P}w_fw_s\mathbf{f}_m(\mathbf{r}_f)\cdot\left[\hat{n}_f\times((\mathbf{r}_f-\mathbf{r}_s)/R^3(1+ikR)e^{-ikR}\times\mathbf{f}_n(\mathbf{r}_s))\right]$
					// 循环遍历场三角形面元上的积分点
					for (int numI = 0; numI < gp.N_points_; ++numI) {
						// 计算场RWG基函数三角形面元上的积分点坐标
						fieldPoint[0] = gp.l1[numI] * vertexF1[0] + gp.l2[numI] * vertexF2[0] + gp.l3[numI] * vertexF3[0];
						fieldPoint[1] = gp.l1[numI] * vertexF1[1] + gp.l2[numI] * vertexF2[1] + gp.l3[numI] * vertexF3[1];
						fieldPoint[2] = gp.l1[numI] * vertexF1[2] + gp.l2[numI] * vertexF2[2] + gp.l3[numI] * vertexF3[2];
						// 获取该场积分点权重
						double wF = gp.weight[numI];

						// 循环遍历源三角形面元上的积分点
						for (int numJ = 0; numJ < gp.N_points_; ++numJ) {
							// 计算源RWG基函数三角形面元上的积分点坐标
							sourcePoint[0] = gp.l1[numJ] * vertexS1[0] + gp.l2[numJ] * vertexS2[0] + gp.l3[numJ] * vertexS3[0];
							sourcePoint[1] = gp.l1[numJ] * vertexS1[1] + gp.l2[numJ] * vertexS2[1] + gp.l3[numJ] * vertexS3[1];
							sourcePoint[2] = gp.l1[numJ] * vertexS1[2] + gp.l2[numJ] * vertexS2[2] + gp.l3[numJ] * vertexS3[2];
							// 获取该源积分点权重
							double wS = gp.weight[numJ];
							// 计算选定场RWG基函数三角形面元上的积分点到源RWG基函数三角形面元上的积分点的距离
							// 高斯积分共享边会选到同一点
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
					// 所用公式: $\sum_{f}\sum_{f}4w_{s}w_{s}A_{f}A_{s}\left[\nabla G(\mathbf{r}_{f},\mathbf{r}_{s})\cdot\left((\mathbf{r}_{f}-\mathbf{r}_{s})\times\mathbf{f}_{n}(\mathbf{r}_{s})\cdot(\mathbf{f}_{m}(\mathbf{r}_{f})\times\hat{n}_{f})\right)\right]$
					for (int numI = 0; numI < gp.N_points_; ++numI) {
						// 计算场RWG基函数三角形面元上的积分点坐标
						fieldPoint[0] = gp.l1[numI] * vertexF1[0] + gp.l2[numI] * vertexF2[0] + gp.l3[numI] * vertexF3[0];
						fieldPoint[1] = gp.l1[numI] * vertexF1[1] + gp.l2[numI] * vertexF2[1] + gp.l3[numI] * vertexF3[1];
						fieldPoint[2] = gp.l1[numI] * vertexF1[2] + gp.l2[numI] * vertexF2[2] + gp.l3[numI] * vertexF3[2];
						// 获取该场积分点权重
						double wF = gp.weight[numI];

						// 循环遍历源三角形面元上的积分点
						for (int numJ = 0; numJ < gp.N_points_; ++numJ) {
							// 计算源RWG基函数三角形面元上的积分点坐标
							sourcePoint[0] = gp.l1[numJ] * vertexS1[0] + gp.l2[numJ] * vertexS2[0] + gp.l3[numJ] * vertexS3[0];
							sourcePoint[1] = gp.l1[numJ] * vertexS1[1] + gp.l2[numJ] * vertexS2[1] + gp.l3[numJ] * vertexS3[1];
							sourcePoint[2] = gp.l1[numJ] * vertexS1[2] + gp.l2[numJ] * vertexS2[2] + gp.l3[numJ] * vertexS3[2];
							// 获取该源积分点权重
							double wS = gp.weight[numJ];

							// 计算选定场RWG基函数三角形面元上的积分点到源RWG基函数三角形面元上的积分点的距离
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
						// 计算场RWG基函数三角形面元上的积分点坐标
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

				// 三元运算符，如果是正三角形，则为1，否则为-1
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
			// 计算场RWG基函数三角形面元上的积分点坐标
			fieldPoint[0] = gp.l1[numF] * vertexF1[0] + gp.l2[numF] * vertexF2[0] + gp.l3[numF] * vertexF3[0];
			fieldPoint[1] = gp.l1[numF] * vertexF1[1] + gp.l2[numF] * vertexF2[1] + gp.l3[numF] * vertexF3[1];
			fieldPoint[2] = gp.l1[numF] * vertexF1[2] + gp.l2[numF] * vertexF2[2] + gp.l3[numF] * vertexF3[2];
			// 获取该积分点权重
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