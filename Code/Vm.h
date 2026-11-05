// Vm.h
#pragma once
#ifndef VM_H
#define VM_H

#include "GaussPoints.h"
#include "OverloadAlgo.h"

#include <array>
#include <cstddef>

namespace RHS {
	namespace detail {
		inline std::complex<double> integrateE(const RWGBase& rwg, const gaussPoints& gp,
			const std::complex<double> k_, const double kInc[3], const double eInc[3])
		{
			const double lm = rwg.edgeLength;
			std::complex<double> value(0.0, 0.0);
			std::array<double, 3> r, rho;

			for (int i = 0; i < gp.N_points_; ++i) {
				r = {
					gp.l1[i] * rwg.positiveFace->vertex1->x + gp.l2[i] * rwg.positiveFace->vertex2->x + gp.l3[i] * rwg.positiveFace->vertex3->x,
					gp.l1[i] * rwg.positiveFace->vertex1->y + gp.l2[i] * rwg.positiveFace->vertex2->y + gp.l3[i] * rwg.positiveFace->vertex3->y,
					gp.l1[i] * rwg.positiveFace->vertex1->z + gp.l2[i] * rwg.positiveFace->vertex2->z + gp.l3[i] * rwg.positiveFace->vertex3->z
				};

				rho = {
					r[0] - rwg.freeVertexPositive->x,
					r[1] - rwg.freeVertexPositive->y,
					r[2] - rwg.freeVertexPositive->z
				};
				std::complex<double> phase = k_ * (kInc[0] * r[0] + kInc[1] * r[1] + kInc[2] * r[2]);
				std::complex<double> e_ = exp(-J * phase);
				value += lm * gp.weight[i] * dot(eInc, rho) * e_;

				r = {
					gp.l1[i] * rwg.negativeFace->vertex1->x + gp.l2[i] * rwg.negativeFace->vertex2->x + gp.l3[i] * rwg.negativeFace->vertex3->x,
					gp.l1[i] * rwg.negativeFace->vertex1->y + gp.l2[i] * rwg.negativeFace->vertex2->y + gp.l3[i] * rwg.negativeFace->vertex3->y,
					gp.l1[i] * rwg.negativeFace->vertex1->z + gp.l2[i] * rwg.negativeFace->vertex2->z + gp.l3[i] * rwg.negativeFace->vertex3->z
				};

				rho = {
					rwg.freeVertexNegative->x - r[0],
					rwg.freeVertexNegative->y - r[1],
					rwg.freeVertexNegative->z - r[2]
				};
				phase = k_ * (kInc[0] * r[0] + kInc[1] * r[1] + kInc[2] * r[2]);
				e_ = exp(-J * phase);
				value += lm * gp.weight[i] * dot(eInc, rho) * e_;
			}

			return value;
		}

		inline std::complex<double> integrateM(const RWGBase& rwg, const gaussPoints& gp,
			const std::complex<double> k_, const double kInc[3], const double hInc[3])
		{
			const double lm = rwg.edgeLength;
			std::complex<double> value(0.0, 0.0);
			std::array<double, 3> r, rho;

			const std::array<double, 3> npos = {
				rwg.positiveFace->externalNormalVector.x,
				rwg.positiveFace->externalNormalVector.y,
				rwg.positiveFace->externalNormalVector.z
			};

			const std::array<double, 3> nneg = {
				rwg.negativeFace->externalNormalVector.x,
				rwg.negativeFace->externalNormalVector.y,
				rwg.negativeFace->externalNormalVector.z
			};

			for (int i = 0; i < gp.N_points_; ++i) {
				r = {
					gp.l1[i] * rwg.positiveFace->vertex1->x + gp.l2[i] * rwg.positiveFace->vertex2->x + gp.l3[i] * rwg.positiveFace->vertex3->x,
					gp.l1[i] * rwg.positiveFace->vertex1->y + gp.l2[i] * rwg.positiveFace->vertex2->y + gp.l3[i] * rwg.positiveFace->vertex3->y,
					gp.l1[i] * rwg.positiveFace->vertex1->z + gp.l2[i] * rwg.positiveFace->vertex2->z + gp.l3[i] * rwg.positiveFace->vertex3->z
				};

				rho = {
					r[0] - rwg.freeVertexPositive->x,
					r[1] - rwg.freeVertexPositive->y,
					r[2] - rwg.freeVertexPositive->z
				};
				std::complex<double> phase = k_ * (kInc[0] * r[0] + kInc[1] * r[1] + kInc[2] * r[2]);
				std::complex<double> e_ = exp(-J * phase);
				value += lm * gp.weight[i] * my_math::cross_dot(npos, hInc, rho) * e_;

				r = {
					gp.l1[i] * rwg.negativeFace->vertex1->x + gp.l2[i] * rwg.negativeFace->vertex2->x + gp.l3[i] * rwg.negativeFace->vertex3->x,
					gp.l1[i] * rwg.negativeFace->vertex1->y + gp.l2[i] * rwg.negativeFace->vertex2->y + gp.l3[i] * rwg.negativeFace->vertex3->y,
					gp.l1[i] * rwg.negativeFace->vertex1->z + gp.l2[i] * rwg.negativeFace->vertex2->z + gp.l3[i] * rwg.negativeFace->vertex3->z
				};

				rho = {
					rwg.freeVertexNegative->x - r[0],
					rwg.freeVertexNegative->y - r[1],
					rwg.freeVertexNegative->z - r[2]
				};
				phase = k_ * (kInc[0] * r[0] + kInc[1] * r[1] + kInc[2] * r[2]);
				e_ = exp(-J * phase);
				value += lm * gp.weight[i] * my_math::cross_dot(nneg, hInc, rho) * e_;
			}

			return value;
		}
	} // namespace detail

	inline std::complex<double> computeV_E_entry(const std::vector<RWGBase>& rwgs, const gaussPoints& gp,
		const std::complex<double> k_, const double kInc[3], const double eInc[3], int rwgid)
	{
		return detail::integrateE(rwgs.at(static_cast<std::size_t>(rwgid)), gp, k_, kInc, eInc);
	}

	inline std::complex<double> computeV_M_entry(const std::vector<RWGBase>& rwgs, const gaussPoints& gp,
		const std::complex<double> k_, const double kInc[3], const double hInc[3], int rwgid)
	{
		return detail::integrateM(rwgs.at(static_cast<std::size_t>(rwgid)), gp, k_, kInc, hInc);
	}

	inline std::complex<double> computeV_C_entry(const std::vector<RWGBase>& rwgs, const gaussPoints& gp,
		const std::complex<double> k_, const double kInc[3], const double eInc[3],
		const double hInc[3], const double alpha_, int rwgid)
	{
		const auto& rwg = rwgs.at(static_cast<std::size_t>(rwgid));
		return alpha_ * detail::integrateE(rwg, gp, k_, kInc, eInc) +
			(1.0 - alpha_) * detail::integrateM(rwg, gp, k_, kInc, hInc);
	}

	inline void computeV_E(const std::vector<RWGBase>& rwgs, const gaussPoints& gp,
		const std::complex<double> k_, double kInc[3], double eInc[3],
		std::vector<std::complex<double>>& Vm)
	{
		Vm.resize(rwgs.size());
		for (std::size_t i = 0; i < rwgs.size(); ++i) {
			Vm[i] = detail::integrateE(rwgs[i], gp, k_, kInc, eInc);
		}
	}

	inline void computeV_M(const std::vector<RWGBase>& rwgs, const gaussPoints& gp,
		const std::complex<double> k_, double kInc[3], double hInc[3],
		std::vector<std::complex<double>>& Vm)
	{
		Vm.resize(rwgs.size());
		for (std::size_t i = 0; i < rwgs.size(); ++i) {
			Vm[i] = detail::integrateM(rwgs[i], gp, k_, kInc, hInc);
		}
	}

	inline void computeV_C(const std::vector<RWGBase>& rwgs, const gaussPoints& gausspoint,
		const std::complex<double> k_, double kInc[3], double eInc[3], double hInc[3],
		const double alpha_, std::vector<std::complex<double>>& Vm)
	{
		Vm.resize(rwgs.size());
		for (std::size_t i = 0; i < rwgs.size(); ++i) {
			Vm[i] = alpha_ * detail::integrateE(rwgs[i], gausspoint, k_, kInc, eInc) +
				(1.0 - alpha_) * detail::integrateM(rwgs[i], gausspoint, k_, kInc, hInc);
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
