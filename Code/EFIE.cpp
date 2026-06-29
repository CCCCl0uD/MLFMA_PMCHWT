// EFIE.cpp
#include "EFIE.h"
#include "OverloadAlgo.h"

namespace EFIE {
	void computeEFIE_Zij(const RWGBase* rwgField, const RWGBase* rwgSource, const std::complex<double> k_, const std::complex<double> eta_,
		const gaussPoints& gp, std::complex<double>& Zij_near_E) {
		int rwgFieldTriId, rwgSourceTriId;
		double areaTriF, areaTriS;
		double constRWGField, constRWGSource;
		double* nonPublicEdge_vertexF, * nonPublicEdge_vertexS;
		double* vertexF1, * vertexF2, * vertexF3;
		double* vertexS1, * vertexS2, * vertexS3;
		double* normalVectorS;
		double lF, lS;

		std::vector<double> fieldPoint(3);
		std::vector<double> sourcePoint(3);
		std::vector<double> pomocvF(3);
		std::vector<double> pomocvS(3);
		std::vector<double> RWGf(3);
		std::vector<double> RWGs(3);

		for (int i = 0; i < 2; ++i) {
			lF = rwgField->edgeLength;
			// 获得场RWG基函数的正负三角形面元信息
			if (i == 0) {
				// 如果是负三角形，则获取负三角形的id
				rwgFieldTriId = rwgField->negativeFace->faceId;
				// 场RWG基函数负三角形面积
				areaTriF = rwgField->negativeFace->triangularArea;
				constRWGField = 1.0 / (2.0 * areaTriF);
				// 获取场RWG基函数的非公共边顶点
				// 注意：这里的非公共边顶点是负三角形的顶点
				nonPublicEdge_vertexF = &(rwgField->freeVertexNegative->x);
				// 获得场RWG基函数的三顶点坐标
				vertexF1 = &(rwgField->negativeFace->vertex1->x);
				vertexF2 = &(rwgField->negativeFace->vertex2->x);
				vertexF3 = &(rwgField->negativeFace->vertex3->x);
			}
			else
			{
				// 如果是正三角形，则获取正三角形的id
				rwgFieldTriId = rwgField->positiveFace->faceId;
				// 场RWG基函数正三角形面积
				areaTriF = rwgField->positiveFace->triangularArea;
				constRWGField = 1.0 / (2.0 * areaTriF);
				// 获取场RWG基函数的非公共边顶点
				nonPublicEdge_vertexF = &(rwgField->freeVertexPositive->x);
				// 获得场RWG基函数的三顶点坐标
				vertexF1 = &(rwgField->positiveFace->vertex1->x);
				vertexF2 = &(rwgField->positiveFace->vertex2->x);
				vertexF3 = &(rwgField->positiveFace->vertex3->x);
			}
			for (int j = 0; j < 2; ++j) {
				lS = rwgSource->edgeLength;
				// 获得源RWG基函数的正负三角形面元信息
				if (j == 0) {
					// 如果是负三角形，则获取负三角形的id
					rwgSourceTriId = rwgSource->negativeFace->faceId;
					// 源RWG基函数负三角形面积
					areaTriS = rwgSource->negativeFace->triangularArea;
					constRWGSource = 1.0 / (2.0 * areaTriS);
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
				else
				{
					// 如果是正三角形，则获取正三角形的id
					rwgSourceTriId = rwgSource->positiveFace->faceId;
					//  源RWG基函数正三角形面积
					areaTriS = rwgSource->positiveFace->triangularArea;
					constRWGSource = 1.0 / (2.0 * areaTriS);
					// 获取源RWG基函数的非公共边顶点
					nonPublicEdge_vertexS = &(rwgSource->freeVertexPositive->x);
					// 获得源RWG基函数的三顶点坐标
					vertexS1 = &(rwgSource->positiveFace->vertex1->x);
					vertexS2 = &(rwgSource->positiveFace->vertex2->x);
					vertexS3 = &(rwgSource->positiveFace->vertex3->x);
					// 获得源RWG基函数正三角形的法向量
					normalVectorS = &(rwgSource->positiveFace->externalNormalVector.x);
				}
				// 初始化结果
				std::complex<double> alok(0.0, 0.0);
				std::complex<double> alok1(0.0, 0.0);

				if (rwgFieldTriId != rwgSourceTriId) {
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

							// 计算格林函数的值
							std::complex<double> G = exp(-J * k_ * R) / R;

							// 计算场RWG基函数三角形面元积分点到该面元非公共棱边顶点的距离
							sub(&pomocvF[0], &fieldPoint[0], nonPublicEdge_vertexF);

							// 计算该场面元下的RWG基函数值
							mul(&RWGf[0], &pomocvF[0], constRWGField);

							// 计算源RWG基函数三角形面元积分点到该面元非公共棱边顶点的距离
							sub(&pomocvS[0], &sourcePoint[0], nonPublicEdge_vertexS);

							// 计算该源面元下的RWG基函数值
							mul(&RWGs[0], &pomocvS[0], constRWGSource);

							alok1 += wF * wS * 4.0 * areaTriF * areaTriS *
								(J * k_ * eta_ * G * dot(&RWGf[0], &RWGs[0]) -
									((J * eta_) / k_) * G * (1.0 / (areaTriF * areaTriS)));
						}
					}
				}
				else {
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
							double R = norm(&fieldPoint[0], &sourcePoint[0]);
							// 计算格林函数的值
							std::complex<double> G_Sing;
							if (R < eps0) {
								G_Sing = -J * k_;
							}
							else {
								G_Sing = ((exp(-J * k_ * R) / R) - 1.0 / R);
							}
							// 计算场RWG基函数三角形面元积分点到该面元非公共棱边顶点的距离
							sub(&pomocvF[0], &fieldPoint[0], nonPublicEdge_vertexF);
							// 场RWG基函数常量
							double constRWGField = (1.0 / (2.0 * areaTriF));
							// 计算该场面元下的RWG基函数值
							mul(&RWGf[0], &pomocvF[0], constRWGField);

							// 计算源RWG基函数三角形面元积分点到该面元非公共棱边顶点的距离
							sub(&pomocvS[0], &sourcePoint[0], nonPublicEdge_vertexS);
							// 源RWG基函数常量
							double constRWGSource = (1.0 / (2.0 * areaTriS));
							// 计算该源面元下的RWG基函数值
							mul(&RWGs[0], &pomocvS[0], constRWGSource);

							alok1 += wF * wS * 4.0 * areaTriF * areaTriS *
								(J * k_ * eta_ * G_Sing * dot(&RWGf[0], &RWGs[0]) -
									((J * eta_) / k_) * G_Sing * (1.0 / (areaTriF * areaTriS)));
						}
					}
					singularityEFIE(areaTriF, areaTriS,
						vertexF1, vertexF2, vertexF3, nonPublicEdge_vertexF, constRWGField,
						vertexS1, vertexS2, vertexS3, nonPublicEdge_vertexS, constRWGSource, normalVectorS,
						gp, k_, eta_, alok);
					alok1 += alok;
				}
				// 三元运算符，如果是正三角形，则为1，否则为-1
				const int signSymbolF = (i == 0) ? -1 : 1;
				const int signSymbolS = (j == 0) ? -1 : 1;
				Zij_near_E += (lF * lS) / (4.0 * Pi) * double(signSymbolF * signSymbolS) * alok1;
			}
		}
	}

	void singularityEFIE(double& areaTriF, double& areaTriS,
		double* vertexF1, double* vertexF2, double* vertexF3, double* nonPublicEdge_vertexF, const double recAf,
		double* vertexS1, double* vertexS2, double* vertexS3, double* nonPublicEdge_vertexS, const double recAs, double* normalVectorS,
		const gaussPoints& gp, const std::complex<double> k_, const std::complex<double> eta_, std::complex<double>& alok)
	{
		std::vector<double> fieldPoint(3);

		std::vector<double> npom1(3);
		std::vector<double> npom2(3);
		std::vector<double> npom3(3);
		std::vector<double> rho(3);
		std::vector<double> rho0N(3);
		std::vector<double> rhotmp(3);
		std::vector<double> rho1(3);
		std::vector<double> rho2(3);
		std::vector<double> rho3(3);
		std::vector<double> KK(3);
		std::vector<double> ll1(3);
		std::vector<double> ll2(3);
		std::vector<double> ll3(3);
		std::vector<double> u1(3);
		std::vector<double> u2(3);
		std::vector<double> u3(3);
		std::vector<double> vecpom(3);
		std::vector<double> p1nulvec(3);
		std::vector<double> p2nulvec(3);
		std::vector<double> p3nulvec(3);
		std::vector<double> pomVEC(3);
		std::vector<double> pomocvF(3);
		std::vector<double> pomocvS(3);
		std::vector<double> pomNvec(3);
		std::vector<double> pomVec1(3);
		std::vector<double> pomVec2(3);
		std::vector<double> Lpl(3);
		std::vector<double> Lmn(3);
		std::vector<double> Pnul(3);
		std::vector<double> Ppl(3);
		std::vector<double> Pmn(3);
		std::vector<double> D(3);
		std::vector<double> Rnul(3);
		std::vector<double> Rpl(3);
		std::vector<double> Rmn(3);
		std::vector<double> POM(3);
		std::vector<double> POM1(3);
		std::vector<double> POM2(3);
		std::vector<double> Ivec(3);
		std::vector<double> RWGf(3);
		std::vector<double> RWGs(3);
		std::vector<double> Ikon(3);

		std::vector<double*> U(3);
		std::vector<double*> Pnulvec(3);

		double Isca;
		double CONST1;
		double CONST2;
		double LOGpl;
		double LOGmn;
		double ATANpl;
		double ATANmn;
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
			double R1mn = sqrt(pow(p1mn, 2) + pow(d1, 2));
			double R2pl = sqrt(pow(p2pl, 2) + pow(d2, 2));
			double R2mn = sqrt(pow(p2mn, 2) + pow(d2, 2));
			double R3pl = sqrt(pow(p3pl, 2) + pow(d3, 2));
			double R3mn = sqrt(pow(p3mn, 2) + pow(d3, 2));

			Rpl[0] = R1pl;
			Rpl[1] = R2pl;
			Rpl[2] = R3pl;

			Rmn[0] = R1mn;
			Rmn[1] = R2mn;
			Rmn[2] = R3mn;

			Ivec[0] = 0.0;
			Ivec[1] = 0.0;
			Ivec[2] = 0.0;

			Isca = 0.0;// 奇异积分标量部分

			/*计算奇异积分*/
			// 矢量部分
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

				CONST1 = 0.5 * (pow(Rnul[iv], 2) * (LOGpl - LOGmn) + Rpl[iv] * Lpl[iv] - Rmn[iv] * Lmn[iv]);

				mul(&pomVEC[0], U[iv], CONST1);
				add(&Ivec[0], &Ivec[0], &pomVEC[0]);
			}

			// 标量部分
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
				Isca = Isca + dot(Pnulvec[is], U[is]) * (Pnul[is] * (LOGpl - LOGmn) - abs(D[is]) * (ATANpl - ATANmn));
			}

			// 伽辽金测试
			sub(&pomocvF[0], &fieldPoint[0], nonPublicEdge_vertexF);
			mul(&RWGf[0], &pomocvF[0], recAf);

			double K1 = dot(&normalVectorS[0], nonPublicEdge_vertexS);
			mul(&pomNvec[0], &normalVectorS[0], K1);
			sub(&rho0N[0], nonPublicEdge_vertexS, &pomNvec[0]);

			sub(&rhotmp[0], &rho[0], &rho0N[0]);
			mul(&pomVec1[0], &rhotmp[0], Isca);
			add(&pomVec2[0], &Ivec[0], &pomVec1[0]);
			mul(&Ikon[0], &pomVec2[0], recAs);
			alok += wF * 2.0 * areaTriF * (J * k_ * eta_ * dot(&RWGf[0], &Ikon[0])
				- ((J * eta_) / k_) * (double(1.0 / areaTriF)) * (double(1.0 / areaTriS)) * (Isca));
		}
	}
}