// Vm.h
#pragma once
#ifndef VM_H
#define VM_H

#include "GaussPoints.h"
#include "OverloadAlgo.h"

namespace RHS {
	inline void computeV_E(const std::vector<RWGBase>& rwgs, const gaussPoints& gp,
		const std::complex<double> k_, double kInc[3], double eInc[3],
		std::vector<std::complex<double>>& Vm)
	{
		int row = static_cast<int>(rwgs.size());
		Vm.resize(row);
		fill(Vm.begin(), Vm.end(), std::complex<double>(0.0));

		for (int rwgid = 0; rwgid < row; ++rwgid) {
			const RWGBase* rwgTemp = &rwgs[rwgid];
			double lm = rwgTemp->edgeLength;

			std::array<double, 3> npos = {
				rwgTemp->positiveFace->externalNormalVector.x,
				rwgTemp->positiveFace->externalNormalVector.y,
				rwgTemp->positiveFace->externalNormalVector.z
			};

			std::array<double, 3> nneg = {
				rwgTemp->negativeFace->externalNormalVector.x,
				rwgTemp->negativeFace->externalNormalVector.y,
				rwgTemp->negativeFace->externalNormalVector.z
			};

			std::array<double, 3> r, rho;

			for (int i = 0; i < gp.N_points_; ++i) {
				r = {
					gp.l1[i] * rwgTemp->positiveFace->vertex1->x + gp.l2[i] * rwgTemp->positiveFace->vertex2->x + gp.l3[i] * rwgTemp->positiveFace->vertex3->x,
					gp.l1[i] * rwgTemp->positiveFace->vertex1->y + gp.l2[i] * rwgTemp->positiveFace->vertex2->y + gp.l3[i] * rwgTemp->positiveFace->vertex3->y,
					gp.l1[i] * rwgTemp->positiveFace->vertex1->z + gp.l2[i] * rwgTemp->positiveFace->vertex2->z + gp.l3[i] * rwgTemp->positiveFace->vertex3->z
				};

				rho = {
					r[0] - rwgTemp->freeVertexPositive->x,
					r[1] - rwgTemp->freeVertexPositive->y,
					r[2] - rwgTemp->freeVertexPositive->z
				};
				std::complex<double> phase = k_ * (kInc[0] * r[0] + kInc[1] * r[1] + kInc[2] * r[2]);
				std::complex<double> e_ = exp(-J * phase);

				Vm[rwgid] = Vm[rwgid] + lm * gp.weight[i] * dot(eInc, rho) * e_;

				r = {
					gp.l1[i] * rwgTemp->negativeFace->vertex1->x + gp.l2[i] * rwgTemp->negativeFace->vertex2->x + gp.l3[i] * rwgTemp->negativeFace->vertex3->x,
					gp.l1[i] * rwgTemp->negativeFace->vertex1->y + gp.l2[i] * rwgTemp->negativeFace->vertex2->y + gp.l3[i] * rwgTemp->negativeFace->vertex3->y,
					gp.l1[i] * rwgTemp->negativeFace->vertex1->z + gp.l2[i] * rwgTemp->negativeFace->vertex2->z + gp.l3[i] * rwgTemp->negativeFace->vertex3->z
				};

				rho = {
					rwgTemp->freeVertexNegative->x - r[0],
					rwgTemp->freeVertexNegative->y - r[1],
					rwgTemp->freeVertexNegative->z - r[2]
				};

				phase = k_ * (kInc[0] * r[0] + kInc[1] * r[1] + kInc[2] * r[2]);
				e_ = exp(-J * phase);

				Vm[rwgid] = Vm[rwgid] + lm * gp.weight[i] * dot(eInc, rho) * e_;
			}
		}
	}

	inline void computeV_M(const std::vector<RWGBase>& rwgs, const gaussPoints& gp,
		const std::complex<double> k_, double kInc[3], double hInc[3],
		std::vector<std::complex<double>>& Vm)
	{
		int row = static_cast<int>(rwgs.size());
		Vm.assign(row, std::complex<double>(0.0));

		for (int rwgid = 0; rwgid < row; ++rwgid) {
			const RWGBase* rwgTemp = &rwgs[rwgid];

			double lm = rwgTemp->edgeLength;

			std::array<double, 3> npos = {
				rwgTemp->positiveFace->externalNormalVector.x,
				rwgTemp->positiveFace->externalNormalVector.y,
				rwgTemp->positiveFace->externalNormalVector.z
			};

			std::array<double, 3> nneg = {
				rwgTemp->negativeFace->externalNormalVector.x,
				rwgTemp->negativeFace->externalNormalVector.y,
				rwgTemp->negativeFace->externalNormalVector.z
			};

			std::array<double, 3> r, rho;

			for (int i = 0; i < gp.N_points_; ++i) {
				r = {
					gp.l1[i] * rwgTemp->positiveFace->vertex1->x + gp.l2[i] * rwgTemp->positiveFace->vertex2->x + gp.l3[i] * rwgTemp->positiveFace->vertex3->x,
					gp.l1[i] * rwgTemp->positiveFace->vertex1->y + gp.l2[i] * rwgTemp->positiveFace->vertex2->y + gp.l3[i] * rwgTemp->positiveFace->vertex3->y,
					gp.l1[i] * rwgTemp->positiveFace->vertex1->z + gp.l2[i] * rwgTemp->positiveFace->vertex2->z + gp.l3[i] * rwgTemp->positiveFace->vertex3->z
				};

				rho = {
					r[0] - rwgTemp->freeVertexPositive->x,
					r[1] - rwgTemp->freeVertexPositive->y,
					r[2] - rwgTemp->freeVertexPositive->z
				};
				std::complex<double> phase = k_ * (kInc[0] * r[0] + kInc[1] * r[1] + kInc[2] * r[2]);
				std::complex<double> e_ = exp(-J * phase);

				Vm[rwgid] = Vm[rwgid] + lm * gp.weight[i] * my_math::cross_dot(npos, hInc, rho) * e_;

				r = {
					gp.l1[i] * rwgTemp->negativeFace->vertex1->x + gp.l2[i] * rwgTemp->negativeFace->vertex2->x + gp.l3[i] * rwgTemp->negativeFace->vertex3->x,
					gp.l1[i] * rwgTemp->negativeFace->vertex1->y + gp.l2[i] * rwgTemp->negativeFace->vertex2->y + gp.l3[i] * rwgTemp->negativeFace->vertex3->y,
					gp.l1[i] * rwgTemp->negativeFace->vertex1->z + gp.l2[i] * rwgTemp->negativeFace->vertex2->z + gp.l3[i] * rwgTemp->negativeFace->vertex3->z
				};

				rho = {
					rwgTemp->freeVertexNegative->x - r[0],
					rwgTemp->freeVertexNegative->y - r[1],
					rwgTemp->freeVertexNegative->z - r[2]
				};

				phase = k_ * (kInc[0] * r[0] + kInc[1] * r[1] + kInc[2] * r[2]);
				e_ = exp(-J * phase);

				Vm[rwgid] = Vm[rwgid] + lm * gp.weight[i] * my_math::cross_dot(nneg, hInc, rho) * e_;
			}
		}
	}

	inline void computeV_C(
		const std::vector<RWGBase>& rwgs, const gaussPoints& gausspoint, const std::complex<double> k_,
		double kInc[3], double eInc[3], double hInc[3], const double alpha_,
		std::vector<std::complex<double>>& Vm) {
		std::vector<std::complex<double>> V1, V2;
		computeV_E(rwgs, gausspoint, k_, kInc, eInc, V1);
		computeV_M(rwgs, gausspoint, k_, kInc, hInc, V2);
		int row = static_cast<int>(rwgs.size());
		Vm.resize(row);
		for (int i = 0; i < row; ++i) {
			Vm[i] = alpha_ * V1[i] + (1.0 - alpha_) * V2[i];
		}
	}

	inline void computeV_PMCHWT(
		const std::vector<RWGBase>& rwgs, const gaussPoints& gp,
		const std::complex<double> k1,
		double kInc[3], double eInc[3], double hInc[3],
		std::vector<std::complex<double>>& Vm)
	{
		int N = static_cast<int>(rwgs.size());
		Vm.assign(2 * N, std::complex<double>(0.0));

		for (int rwgid = 0; rwgid < N; ++rwgid) {
			const RWGBase* r = &rwgs[rwgid];
			double lm = r->edgeLength;

			std::array<double, 3> rr, rho;

			for (int i = 0; i < gp.N_points_; ++i) {
				rr = {
					gp.l1[i] * r->positiveFace->vertex1->x + gp.l2[i] * r->positiveFace->vertex2->x + gp.l3[i] * r->positiveFace->vertex3->x,
					gp.l1[i] * r->positiveFace->vertex1->y + gp.l2[i] * r->positiveFace->vertex2->y + gp.l3[i] * r->positiveFace->vertex3->y,
					gp.l1[i] * r->positiveFace->vertex1->z + gp.l2[i] * r->positiveFace->vertex2->z + gp.l3[i] * r->positiveFace->vertex3->z
				};

				rho = {
					rr[0] - r->freeVertexPositive->x,
					rr[1] - r->freeVertexPositive->y,
					rr[2] - r->freeVertexPositive->z
				};

				std::complex<double> phase1 = k1 * (kInc[0] * rr[0] + kInc[1] * rr[1] + kInc[2] * rr[2]);
				std::complex<double> e_k1 = exp(-J * phase1);

				Vm[rwgid] += lm * gp.weight[i] * dot(eInc, rho) * e_k1;
				Vm[N + rwgid] += lm * gp.weight[i] * dot(hInc, rho) * e_k1;

				rr = {
					gp.l1[i] * r->negativeFace->vertex1->x + gp.l2[i] * r->negativeFace->vertex2->x + gp.l3[i] * r->negativeFace->vertex3->x,
					gp.l1[i] * r->negativeFace->vertex1->y + gp.l2[i] * r->negativeFace->vertex2->y + gp.l3[i] * r->negativeFace->vertex3->y,
					gp.l1[i] * r->negativeFace->vertex1->z + gp.l2[i] * r->negativeFace->vertex2->z + gp.l3[i] * r->negativeFace->vertex3->z
				};

				rho = {
					r->freeVertexNegative->x - rr[0],
					r->freeVertexNegative->y - rr[1],
					r->freeVertexNegative->z - rr[2]
				};

				phase1 = k1 * (kInc[0] * rr[0] + kInc[1] * rr[1] + kInc[2] * rr[2]);
				e_k1 = exp(-J * phase1);

				Vm[rwgid] += lm * gp.weight[i] * dot(eInc, rho) * e_k1;
				Vm[N + rwgid] += lm * gp.weight[i] * dot(hInc, rho) * e_k1;
			}
		}
	}
}

#endif // !VM_H
