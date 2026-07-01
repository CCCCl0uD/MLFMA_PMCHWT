// FMM.cpp
#include "FMM.h"
#include "AlgoFMM.h"
#include "GMRES.h"
#include "CGS.h"
#include "ProcessFMM.h"
#include "Zmartix.h"
#include "Vm.h"
#include "RCS.h"

FMM::FMM(const RCSExportConfig& cfg, const int selectIntegralEqu, const int selectMatrixSolver, const std::string selectMono_Dual, const std::string pol_wave,
	const std::vector<std::vector<OCTree::Node*>>& octreeNodes,
	const std::vector<RWGBase>& rwgs, const int maxLevel_,
	const gaussPoints& gausspoint, const double E0, const EMSource& wave)
	: octreeNodes_(octreeNodes), rwgs(rwgs), maxLevel_(maxLevel_), gausspoint(gausspoint), E0(E0), wave(wave), integralEquType_(selectIntegralEqu), matrixSolverType_(selectMatrixSolver)
{
	// Compute base data: k, top level theta, L, levelSpan
	// Compute Vsmi/Vfmj matrices theta/phi angles and weights from maxLevel_ down to level 2
	AlgoFMM::computeBase(octreeNodes_, wave, integralEquType_, maxLevel_, static_cast<int>(rwgs.size()), L_k1, L_k2, row, levelSpan,
		WGL_k1, WGL_k2, WGL_phi_k1, WGL_phi_k2, theta_level_k1, theta_level_k2, phi_level_k1, phi_level_k2, kp_lvl_k1, kp_lvl_k2);

	if (integralEquType_ != 2) {
		AlgoFMM::computeFMMFarGroups(octreeNodes_[maxLevel_], farVec, node_far, node_far_id, node_farVec_id);

		AlgoFMM::computeFMM_far(farVec, node_far, node_far_id, node_farVec_id,
			TF_fmm, kp_lvl_k1, wave.k1(), L_k1);

		if (integralEquType_ == 0) {
			// k1: Vsmi/Vfmj
			ProcessFMM::computeVsmi_Vfmj_Efie(Vsmi, Vfmj, octreeNodes_[maxLevel_], kp_lvl_k1[0], wave,
				row, gausspoint, alpha, 1);// Compute Vsmi and Vfmj (EFIE)

			// PMCHWT: Z_near
			Zmartix::computeZ_fmm_pec_efie(octreeNodes_[maxLevel_], rwgs, Z_near, Z_near_id, row, wave,
				gausspoint, 1);

			if (selectMono_Dual == "dual") {
				fmm_Dual_Pec_Efie(cfg, pol_wave);
			}

			if (selectMono_Dual == "mono") {
				fmm_Mono_Pec_Efie(cfg, pol_wave);
			}
		}

		if (integralEquType_ == 1) {
			// k1: Vsmi/Vfmj
			ProcessFMM::computeVsmi_Vfmj_Cfie(Vsmi, Vfmj, octreeNodes_[maxLevel_], kp_lvl_k1[0], wave,
				row, gausspoint, alpha, 1);

			// PMCHWT: Z_near
			Zmartix::computeZ_fmm_pec_cfie(octreeNodes_[maxLevel_], rwgs, Z_near, Z_near_id, row, wave,
				gausspoint, alpha, 1);

			if (selectMono_Dual == "dual") {
				fmm_Dual_Pec_Cfie(cfg, pol_wave);
			}

			if (selectMono_Dual == "mono") {
				fmm_Mono_Pec_Cfie(cfg, pol_wave);
			}
		}
	}
	else {
		AlgoFMM::computeFMMFarGroups(octreeNodes_[maxLevel_], farVec, node_far, node_far_id, node_farVec_id);

		// k1: TF_fmm1
		AlgoFMM::computeFMM_far(farVec, node_far, node_far_id, node_farVec_id,
			TF_fmm, kp_lvl_k1, wave.k1(), L_k1);

		// k2: TF_fmm2
		AlgoFMM::computeFMM_far(farVec, node_far, node_far_id, node_farVec_id,
			TF_fmm2, kp_lvl_k2, wave.k2(), L_k2);

		// k1: Vsmi1/Vfmj1
		ProcessFMM::computeVsmi_Vfmj_Die(Vsmi, Vfmj,
			octreeNodes_[maxLevel_], kp_lvl_k1[0], wave.k1(), wave.eta1(),
			row, gausspoint, alpha, 1);

		// k2: Vsmi2/Vfmj2
		ProcessFMM::computeVsmi_Vfmj_Die(Vsmi2, Vfmj2,
			octreeNodes_[maxLevel_], kp_lvl_k2[0], wave.k2(), wave.eta2(),
			row, gausspoint, alpha, 1);

		// PMCHWT: Z_near
		Zmartix::computeZ_fmm_die_pmchwt(octreeNodes_[maxLevel_], rwgs,
			Z_near, Z_near_id, row, wave, gausspoint, 1);

		if (selectMono_Dual == "dual") {
			fmm_Dual_Die_Pmchwt(cfg, pol_wave);
		}
		if (selectMono_Dual == "mono") {
			fmm_Mono_Die_Pmchwt(cfg, pol_wave);
		}
	}
}

size_t FMM::computeMem() {
	size_t memBase = memory2D<kp_Point>(kp_lvl_k1) + memory2D<double>(theta_level_k1) + memory2D<double>(phi_level_k1) +
		memory2D<double>(WGL_k1) + memory1D<double>(WGL_phi_k1);
	size_t memVsmi = memory2D<Complex3D>(Vsmi) + memory2D<Complex3D>(Vfmj);
	size_t memZ = memory2D<std::complex<double>>(Z_near) + memory2D<int>(Z_near_id);
	size_t memTF = memory2D<std::complex<double>>(TF_fmm);
	size_t memSmBm = memory2D<Complex3D>(Sm_fmm) + memory2D<Complex3D>(Bm_fmm);
	if (integralEquType_ == 2) {
		memBase += memory2D<kp_Point>(kp_lvl_k2) + memory2D<double>(theta_level_k2) + memory2D<double>(phi_level_k2) +
			memory2D<double>(WGL_k2) + memory1D<double>(WGL_phi_k2);
		memVsmi += memory2D<Complex3D>(Vsmi2) + memory2D<Complex3D>(Vfmj2);
		memTF += memory2D<std::complex<double>>(TF_fmm2);
	}
	return memBase + memVsmi + memZ + memTF + memSmBm;
}

void FMM::matrix_solver(int n, std::complex<double> x[], std::complex<double> rhs[], int itr_max, int mr, double tol_abs, double tol_rel) {
	const int EquType_ = integralEquType_;

	auto ax = [this, EquType_](int n, const std::complex<double>* x, std::complex<double>* w) {
		for (int i = 0; i < n; ++i) {
			w[i] = std::complex<double>(0.0, 0.0);
		}

		int nodenum = static_cast<int>(octreeNodes_[maxLevel_].size());

		if (EquType_ != 2) {
			const int K1 = static_cast<int>(kp_lvl_k1[0].size());

			// ========== PEC ==========
			Sm_fmm.assign(nodenum, std::vector<Complex3D>(K1, Complex3D()));
			Bm_fmm.assign(nodenum, std::vector<Complex3D>(K1, Complex3D()));

#pragma omp parallel
			{
				// Multipole expansion (Mexp)
#pragma omp for schedule(dynamic)
				for (int nodeN = 0; nodeN < nodenum; ++nodeN) {
					OCTree::Node* node_mexp = octreeNodes_[maxLevel_][nodeN];

					const int rwgNum = node_mexp->rwgIndices.size();

					for (int rwg_i = 0; rwg_i < rwgNum; ++rwg_i) {
						int id = node_mexp->rwgIndices[rwg_i]->rwgid;

						for (int kpi = 0; kpi < kp_lvl_k1[0].size(); ++kpi) {
							Sm_fmm[nodeN][kpi].x += Vsmi[id][kpi].x * x[id];
							Sm_fmm[nodeN][kpi].y += Vsmi[id][kpi].y * x[id];
							Sm_fmm[nodeN][kpi].z += Vsmi[id][kpi].z * x[id];
						}
					}
				}
#pragma omp barrier
				// Multipole Expansion to Local Expansion (ML) Bm
#pragma omp for schedule(dynamic)
				for (int nodeN = 0; nodeN < nodenum; ++nodeN) {
					const int farnum = node_far[nodeN].size();

					for (int farn = 0; farn < farnum; ++farn) {
						const auto& TF = TF_fmm[node_farVec_id[nodeN][farn]];
						int nodeid = node_far_id[nodeN][farn];

						for (int kpi = 0; kpi < kp_lvl_k1[0].size(); ++kpi) {
							Bm_fmm[nodeN][kpi].x += TF[kpi] * Sm_fmm[nodeid][kpi].x;
							Bm_fmm[nodeN][kpi].y += TF[kpi] * Sm_fmm[nodeid][kpi].y;
							Bm_fmm[nodeN][kpi].z += TF[kpi] * Sm_fmm[nodeid][kpi].z;
						}
					}
				}
			}
			// Local Expansion (Lexp)
			const int tN = static_cast<int>(theta_level_k1[0].size());
			const int pN = static_cast<int>(phi_level_k1[0].size());
#pragma omp parallel for schedule(dynamic)
			for (int nodeIdx = 0; nodeIdx < nodenum; ++nodeIdx) {
				OCTree::Node* node = octreeNodes_[maxLevel_][nodeIdx];
				const int rwgNum = node->rwgIndices.size();

				for (int rwgCount = 0; rwgCount < rwgNum; ++rwgCount) {
					int rwgid = node->rwgIndices[rwgCount]->rwgid;

					const auto& V = Vfmj[rwgid];
					const auto& B = Bm_fmm[nodeIdx];

					std::complex<double> val(0.0, 0.0);

					for (int tt = 0; tt < tN; ++tt) {
						for (int pp = 0; pp < pN; ++pp) {
							int idx = tt * pN + pp;
							double d2k = WGL_k1[0][tt] * WGL_phi_k1[0];
							val += d2k * dot(V[idx], B[idx]);
						}
					}

					val *= ((wave.k1() * wave.eta1()) / (4.0 * Pi)) * (wave.k1() / (4.0 * Pi));

					for (int j = 0; j < static_cast<int>(Z_near[rwgid].size()); ++j) {
						val += x[Z_near_id[rwgid][j]] * Z_near[rwgid][j];
					}

					w[rwgid] = val;
				}
			}
		}
		else {
			// ===== PMCHWT matrix-vector product =====
			// Unknown ordering:
			//   x[0:N)     = J coefficients
			//   x[N:2*N)   = M coefficients
			//
			// Equation ordering:
			//   w[0:N)     = electric-current testing block
			//   w[N:2*N)   = magnetic-current testing block

			const int N = row;
			if (n != 2 * N) {
				throw std::invalid_argument("FMM PMCHWT matrix_solver requires n == 2 * row.");
			}

			const int K1 = static_cast<int>(kp_lvl_k1[0].size());
			const int K2 = static_cast<int>(kp_lvl_k2[0].size());

			const std::complex<double> etai = wave.etai();

			const int tN1 = static_cast<int>(theta_level_k1[0].size());
			const int pN1 = static_cast<int>(phi_level_k1[0].size());
			const int tN2 = static_cast<int>(theta_level_k2[0].size());
			const int pN2 = static_cast<int>(phi_level_k2[0].size());

			// const参数存疑
			const std::complex<double> const1 = wave.k1() * wave.k1() / (16.0 * Pi * Pi);
			const std::complex<double> const2 = wave.k2() * wave.k2() / (16.0 * Pi * Pi);

			auto aggregate = [&](const std::vector<std::vector<Complex3D>>& Vsm,
				int xOffset, int K, std::vector<std::vector<Complex3D>>& Sm) {
#pragma omp parallel for schedule(dynamic)
					for (int nodeIdx = 0; nodeIdx < nodenum; ++nodeIdx) {
						OCTree::Node* node = octreeNodes_[maxLevel_][nodeIdx];
						const int rwgNum = static_cast<int>(node->rwgIndices.size());

						for (int r = 0; r < rwgNum; ++r) {
							const int rwgid = node->rwgIndices[r]->rwgid;
							const std::complex<double> coeff = x[xOffset + rwgid];

							for (int kpi = 0; kpi < K; ++kpi) {
								Sm[nodeIdx][kpi].x += Vsm[rwgid][kpi].x * coeff;
								Sm[nodeIdx][kpi].y += Vsm[rwgid][kpi].y * coeff;
								Sm[nodeIdx][kpi].z += Vsm[rwgid][kpi].z * coeff;
							}
						}
					}
				};

			auto translate = [&](const std::vector<std::vector<std::complex<double>>>& TF_all,
				int K, const std::vector<std::vector<Complex3D>>& Sm,
				std::vector<std::vector<Complex3D>>& Bm) {
#pragma omp parallel for schedule(dynamic)
					for (int nodeIdx = 0; nodeIdx < nodenum; ++nodeIdx) {
						const int farnum = static_cast<int>(node_far[nodeIdx].size());

						for (int f = 0; f < farnum; ++f) {
							const int farId = node_far_id[nodeIdx][f];
							const auto& TF = TF_all[node_farVec_id[nodeIdx][f]];

							for (int kpi = 0; kpi < K; ++kpi) {
								Bm[nodeIdx][kpi].x += TF[kpi] * Sm[farId][kpi].x;
								Bm[nodeIdx][kpi].y += TF[kpi] * Sm[farId][kpi].y;
								Bm[nodeIdx][kpi].z += TF[kpi] * Sm[farId][kpi].z;
							}
						}
					}
				};

			auto integrateL = [](const std::vector<Complex3D>& V,
				const std::vector<Complex3D>& B, const std::vector<double>& WGL,
				double WGL_phi, int tN, int pN) -> std::complex<double> {
					std::complex<double> val(0.0, 0.0);
					for (int ti = 0; ti < tN; ++ti) {
						for (int pj = 0; pj < pN; ++pj) {
							const int idx = ti * pN + pj;
							const double d2k = WGL[ti] * WGL_phi;
							val += d2k * dot(V[idx], B[idx]);
						}
					}
					return val;
				};

			auto integrateK = [](const std::vector<Complex3D>& V,
				const std::vector<Complex3D>& B, const std::vector<kp_Point>& kp,
				const std::vector<double>& WGL, double WGL_phi, int tN, int pN) -> std::complex<double> {
					std::complex<double> val(0.0, 0.0);
					for (int ti = 0; ti < tN; ++ti) {
						for (int pj = 0; pj < pN; ++pj) {
							const int idx = ti * pN + pj;
							const double d2k = WGL[ti] * WGL_phi;

							const double khat[3] = { kp[idx].x, kp[idx].y, kp[idx].z };

							Complex3D kCrossV;
							cross(kCrossV, khat, V[idx]);
							val += d2k * dot(kCrossV, B[idx]);
						}
					}
					return val;
				};

			// ---------- k1 external space ----------
			{
				std::vector<std::vector<Complex3D>> Sm_J1(
					nodenum, std::vector<Complex3D>(K1, Complex3D()));
				std::vector<std::vector<Complex3D>> Bm_J1(
					nodenum, std::vector<Complex3D>(K1, Complex3D()));
				std::vector<std::vector<Complex3D>> Sm_M1(
					nodenum, std::vector<Complex3D>(K1, Complex3D()));
				std::vector<std::vector<Complex3D>> Bm_M1(
					nodenum, std::vector<Complex3D>(K1, Complex3D()));

				aggregate(Vsmi, 0, K1, Sm_J1);
				aggregate(Vsmi, N, K1, Sm_M1);

				translate(TF_fmm, K1, Sm_J1, Bm_J1);
				translate(TF_fmm, K1, Sm_M1, Bm_M1);

#pragma omp parallel for schedule(dynamic)
				for (int nodeIdx = 0; nodeIdx < nodenum; ++nodeIdx) {
					OCTree::Node* node = octreeNodes_[maxLevel_][nodeIdx];
					const int rwgNum = static_cast<int>(node->rwgIndices.size());

					for (int r = 0; r < rwgNum; ++r) {
						const int rwgid = node->rwgIndices[r]->rwgid;
						const auto& V1 = Vfmj[rwgid];

						const auto L1J = integrateL(
							V1, Bm_J1[nodeIdx], WGL_k1[0], WGL_phi_k1[0], tN1, pN1);
						const auto K1J = integrateK(
							V1, Bm_J1[nodeIdx], kp_lvl_k1[0], WGL_k1[0], WGL_phi_k1[0], tN1, pN1);
						const auto K1M = integrateK(
							V1, Bm_M1[nodeIdx], kp_lvl_k1[0], WGL_k1[0], WGL_phi_k1[0], tN1, pN1);
						const auto L1M = integrateL(
							V1, Bm_M1[nodeIdx], WGL_k1[0], WGL_phi_k1[0], tN1, pN1);

						w[rwgid] += const1 * L1J;
						w[rwgid] += -const1 * K1M;
						w[N + rwgid] += const1 * K1J;
						w[N + rwgid] += const1 * L1M;
					}
				}
			}

			// ---------- k2 internal space ----------
			{
				std::vector<std::vector<Complex3D>> Sm_J2(
					nodenum, std::vector<Complex3D>(K2, Complex3D()));
				std::vector<std::vector<Complex3D>> Bm_J2(
					nodenum, std::vector<Complex3D>(K2, Complex3D()));
				std::vector<std::vector<Complex3D>> Sm_M2(
					nodenum, std::vector<Complex3D>(K2, Complex3D()));
				std::vector<std::vector<Complex3D>> Bm_M2(
					nodenum, std::vector<Complex3D>(K2, Complex3D()));

				aggregate(Vsmi2, 0, K2, Sm_J2);
				aggregate(Vsmi2, N, K2, Sm_M2);

				translate(TF_fmm2, K2, Sm_J2, Bm_J2);
				translate(TF_fmm2, K2, Sm_M2, Bm_M2);

#pragma omp parallel for schedule(dynamic)
				for (int nodeIdx = 0; nodeIdx < nodenum; ++nodeIdx) {
					OCTree::Node* node = octreeNodes_[maxLevel_][nodeIdx];
					const int rwgNum = static_cast<int>(node->rwgIndices.size());

					for (int r = 0; r < rwgNum; ++r) {
						const int rwgid = node->rwgIndices[r]->rwgid;
						const auto& V2 = Vfmj2[rwgid];

						const auto L2J = integrateL(
							V2, Bm_J2[nodeIdx], WGL_k2[0], WGL_phi_k2[0], tN2, pN2);
						const auto K2J = integrateK(
							V2, Bm_J2[nodeIdx], kp_lvl_k2[0], WGL_k2[0], WGL_phi_k2[0], tN2, pN2);
						const auto K2M = integrateK(
							V2, Bm_M2[nodeIdx], kp_lvl_k2[0], WGL_k2[0], WGL_phi_k2[0], tN2, pN2);
						const auto L2M = integrateL(
							V2, Bm_M2[nodeIdx], WGL_k2[0], WGL_phi_k2[0], tN2, pN2);

						w[rwgid] += const2 * etai * L2J;
						w[rwgid] += -const2 * K2M;
						w[N + rwgid] += const2 * K2J;
						w[N + rwgid] += const2 * L2M / etai;
					}
				}
			}

#pragma omp parallel for schedule(dynamic)
			for (int nodeIdx = 0; nodeIdx < nodenum; ++nodeIdx) {
				OCTree::Node* nd = octreeNodes_[maxLevel_][nodeIdx];
				for (int r = 0; r < static_cast<int>(nd->rwgIndices.size()); ++r) {
					int rwgid = nd->rwgIndices[r]->rwgid;

					for (int j = 0; j < static_cast<int>(Z_near[rwgid].size()); ++j) {
						w[rwgid] += x[Z_near_id[rwgid][j]] * Z_near[rwgid][j];
					}

					for (int j = 0; j < static_cast<int>(Z_near[N + rwgid].size()); ++j) {
						w[N + rwgid] += x[Z_near_id[N + rwgid][j]] * Z_near[N + rwgid][j];
					}
				}
			}
		}
		};// end ax lambda
	bool converged = false;

	switch (matrixSolverType_) {
	case 0:
		std::cout << "Matrix solver: GMRES (FMM)\n";
		converged = MGMRES::mgmres(
			n, x, rhs, itr_max, mr, tol_abs, tol_rel, ax);
		if (!converged) {
			std::cerr << "GMRES did not converge (FMM).\n";
		}
		break;

	case 1:
		std::cout << "Matrix solver: CGS (FMM)\n";
		converged = MCGS::CGS::solve(
			n, x, rhs, itr_max, tol_rel, ax);
		if (!converged) {
			std::cerr
				<< "CGS did not converge or encountered breakdown (FMM).\n";
		}
		break;

	default:
		throw std::invalid_argument(
			"Invalid matrixSolverType_ in FMM "
			"(must be 0 for GMRES or 1 for CGS)");
	}
}

void FMM::fmm_Dual_Pec_Efie(const RCSExportConfig& cfg, const std::string pol_wave)
{
	auto computeV = [](FMM& solver, double kInc[3], double eInc[3], double hInc[3]) {
		(void)hInc;
		RHS::computeV_E(solver.rwgs, solver.gausspoint,
			solver.wave.k1(), kInc, eInc, solver.Vm);
		};
	RCSUtils::computeDualStatic(*this, cfg, computeV);
}

void FMM::fmm_Mono_Pec_Efie(const RCSExportConfig& cfg, const std::string pol_wave)
{
	auto computeV = [](FMM& solver, double kInc[3], double eInc[3], double hInc[3]) {
		(void)hInc;
		RHS::computeV_E(solver.rwgs, solver.gausspoint,
			solver.wave.k1(), kInc, eInc, solver.Vm);
		};
	RCSUtils::computeMonoStatic(*this, cfg, computeV);
}

void FMM::fmm_Dual_Pec_Cfie(const RCSExportConfig& cfg, const std::string pol_wave)
{
	auto computeV = [](FMM& solver, double kInc[3], double eInc[3], double hInc[3]) {
		RHS::computeV_C(solver.rwgs, solver.gausspoint,
			solver.wave.k1(), kInc, eInc, hInc, alpha, solver.Vm);
		};
	RCSUtils::computeDualStatic(*this, cfg, computeV);
}

void FMM::fmm_Mono_Pec_Cfie(const RCSExportConfig& cfg, const std::string pol_wave)
{
	auto computeV = [this](FMM& solver, double kInc[3], double eInc[3], double hInc[3]) {
		RHS::computeV_C(solver.rwgs, solver.gausspoint,
			solver.wave.k1(), kInc, eInc, hInc, alpha, solver.Vm);
		};
	RCSUtils::computeMonoStatic(*this, cfg, computeV);
}

void FMM::fmm_Dual_Die_Pmchwt(const RCSExportConfig& cfg, const std::string pol_wave)
{
	auto computeV = [](FMM& solver, double kInc[3], double eInc[3], double hInc[3]) {
		RHS::computeV_PMCHWT(solver.rwgs, solver.gausspoint,
			solver.wave.k1(), kInc, eInc, hInc, solver.Vm);
		};
	RCSUtils::computeDualStatic_PMCHWT(*this, cfg, computeV);
}

void FMM::fmm_Mono_Die_Pmchwt(const RCSExportConfig& cfg, const std::string pol_wave)
{
	auto computeV = [](FMM& solver, double kInc[3], double eInc[3], double hInc[3]) {
		RHS::computeV_PMCHWT(solver.rwgs, solver.gausspoint,
			solver.wave.k1(), kInc, eInc, hInc, solver.Vm);
		};
	RCSUtils::computeMonoStatic_PMCHWT(*this, cfg, computeV);
}

void FMM::fmm_fullZmn_Die_Pmchwt() {
	// ===== PMCHWT matrix-vector product =====
			// Unknown ordering:
			//   x[0:N)     = J coefficients
			//   x[N:2*N)   = M coefficients
			//
			// Equation ordering:
			//   w[0:N)     = electric-current testing block
			//   w[N:2*N)   = magnetic-current testing block

	const int N = row;

	std::vector<std::complex<double>> x(2 * N);
	for (int i = 0; i < 2 * N; ++i) {
		x[i] = std::complex<double>(1.0, 0.0);
	}

	int nodenum = static_cast<int>(octreeNodes_[maxLevel_].size());

	std::vector<std::vector<std::complex<double>>> Z_full;
	Z_full.resize(2 * N);
	for (int m = 0; m < 2 * N; ++m) {
		Z_full[m].assign(2 * N, std::complex<double>(0.0, 0.0));
	}

	const int K1 = static_cast<int>(kp_lvl_k1[0].size());
	const int K2 = static_cast<int>(kp_lvl_k2[0].size());

	const int tN1 = static_cast<int>(theta_level_k1[0].size());
	const int pN1 = static_cast<int>(phi_level_k1[0].size());
	const int tN2 = static_cast<int>(theta_level_k2[0].size());
	const int pN2 = static_cast<int>(phi_level_k2[0].size());

	// const参数存疑
	const std::complex<double> const1 = wave.k1() * wave.k1() / (16.0 * Pi * Pi);
	const std::complex<double> const2 = wave.k2() * wave.k2() / (16.0 * Pi * Pi);

	auto aggregate = [&](const std::vector<std::vector<Complex3D>>& Vsm,
		int xOffset, int K, std::vector<std::vector<Complex3D>>& Sm) {
#pragma omp parallel for schedule(dynamic)
			for (int nodeIdx = 0; nodeIdx < nodenum; ++nodeIdx) {
				OCTree::Node* node = octreeNodes_[maxLevel_][nodeIdx];
				const int rwgNum = static_cast<int>(node->rwgIndices.size());

				for (int r = 0; r < rwgNum; ++r) {
					const int rwgid = node->rwgIndices[r]->rwgid;
					const std::complex<double> coeff = x[xOffset + rwgid];

					for (int kpi = 0; kpi < K; ++kpi) {
						Sm[nodeIdx][kpi].x += Vsm[rwgid][kpi].x * coeff;
						Sm[nodeIdx][kpi].y += Vsm[rwgid][kpi].y * coeff;
						Sm[nodeIdx][kpi].z += Vsm[rwgid][kpi].z * coeff;
					}
				}
			}
		};

	auto translate = [&](const std::vector<std::vector<std::complex<double>>>& TF_all,
		int K, const std::vector<std::vector<Complex3D>>& Sm,
		std::vector<std::vector<Complex3D>>& Bm) {
#pragma omp parallel for schedule(dynamic)
			for (int nodeIdx = 0; nodeIdx < nodenum; ++nodeIdx) {
				const int farnum = static_cast<int>(node_far[nodeIdx].size());

				for (int f = 0; f < farnum; ++f) {
					const int farId = node_far_id[nodeIdx][f];
					const auto& TF = TF_all[node_farVec_id[nodeIdx][f]];

					for (int kpi = 0; kpi < K; ++kpi) {
						Bm[nodeIdx][kpi].x += TF[kpi] * Sm[farId][kpi].x;
						Bm[nodeIdx][kpi].y += TF[kpi] * Sm[farId][kpi].y;
						Bm[nodeIdx][kpi].z += TF[kpi] * Sm[farId][kpi].z;
					}
				}
			}
		};

	auto integrateL = [](const std::vector<Complex3D>& V,
		const std::vector<Complex3D>& B, const std::vector<double>& WGL,
		double WGL_phi, int tN, int pN) -> std::complex<double> {
			std::complex<double> val(0.0, 0.0);
			for (int ti = 0; ti < tN; ++ti) {
				for (int pj = 0; pj < pN; ++pj) {
					const int idx = ti * pN + pj;
					const double d2k = WGL[ti] * WGL_phi;
					val += d2k * dot(V[idx], B[idx]);
				}
			}
			return val;
		};

	auto integrateK = [](const std::vector<Complex3D>& V,
		const std::vector<Complex3D>& B, const std::vector<kp_Point>& kp,
		const std::vector<double>& WGL, double WGL_phi, int tN, int pN) -> std::complex<double> {
			std::complex<double> val(0.0, 0.0);
			for (int ti = 0; ti < tN; ++ti) {
				for (int pj = 0; pj < pN; ++pj) {
					const int idx = ti * pN + pj;
					const double d2k = WGL[ti] * WGL_phi;

					const double khat[3] = { kp[idx].x, kp[idx].y, kp[idx].z };

					Complex3D kCrossV;
					cross(kCrossV, khat, V[idx]);
					val += d2k * dot(kCrossV, B[idx]);
				}
			}
			return val;
		};

	// ---------- k1 external space ----------
	{
		std::vector<std::vector<Complex3D>> Sm_J1(
			nodenum, std::vector<Complex3D>(K1, Complex3D()));
		std::vector<std::vector<Complex3D>> Bm_J1(
			nodenum, std::vector<Complex3D>(K1, Complex3D()));
		std::vector<std::vector<Complex3D>> Sm_M1(
			nodenum, std::vector<Complex3D>(K1, Complex3D()));
		std::vector<std::vector<Complex3D>> Bm_M1(
			nodenum, std::vector<Complex3D>(K1, Complex3D()));

		aggregate(Vsmi, 0, K1, Sm_J1);
		aggregate(Vsmi, N, K1, Sm_M1);

		translate(TF_fmm, K1, Sm_J1, Bm_J1);
		translate(TF_fmm, K1, Sm_M1, Bm_M1);

#pragma omp parallel for schedule(dynamic)
		for (int nodeIdx = 0; nodeIdx < nodenum; ++nodeIdx) {
			OCTree::Node* node = octreeNodes_[maxLevel_][nodeIdx];
			const int rwgNum = static_cast<int>(node->rwgIndices.size());

			for (int r = 0; r < rwgNum; ++r) {
				const int rwgid = node->rwgIndices[r]->rwgid;
				const auto& V1 = Vfmj[rwgid];

				const auto L1J = integrateL(
					V1, Bm_J1[nodeIdx], WGL_k1[0], WGL_phi_k1[0], tN1, pN1);
				const auto K1J = integrateK(
					V1, Bm_J1[nodeIdx], kp_lvl_k1[0], WGL_k1[0], WGL_phi_k1[0], tN1, pN1);
				const auto K1M = integrateK(
					V1, Bm_M1[nodeIdx], kp_lvl_k1[0], WGL_k1[0], WGL_phi_k1[0], tN1, pN1);
				const auto L1M = integrateL(
					V1, Bm_M1[nodeIdx], WGL_k1[0], WGL_phi_k1[0], tN1, pN1);

				Z_full[rwgNum][rwgid] += const1 * L1J;
				Z_full[rwgNum][N + rwgid] += -const1 * K1J;
				Z_full[rwgNum][rwgid] += const1 * K1M;
				Z_full[rwgNum][N + rwgid] += const1 * L1M;
			}
		}
	}

	// ---------- k2 internal space ----------
	{
		std::vector<std::vector<Complex3D>> Sm_J2(
			nodenum, std::vector<Complex3D>(K2, Complex3D()));
		std::vector<std::vector<Complex3D>> Bm_J2(
			nodenum, std::vector<Complex3D>(K2, Complex3D()));
		std::vector<std::vector<Complex3D>> Sm_M2(
			nodenum, std::vector<Complex3D>(K2, Complex3D()));
		std::vector<std::vector<Complex3D>> Bm_M2(
			nodenum, std::vector<Complex3D>(K2, Complex3D()));

		aggregate(Vsmi2, 0, K2, Sm_J2);
		aggregate(Vsmi2, N, K2, Sm_M2);

		translate(TF_fmm2, K2, Sm_J2, Bm_J2);
		translate(TF_fmm2, K2, Sm_M2, Bm_M2);

#pragma omp parallel for schedule(dynamic)
		for (int nodeIdx = 0; nodeIdx < nodenum; ++nodeIdx) {
			OCTree::Node* node = octreeNodes_[maxLevel_][nodeIdx];
			const int rwgNum = static_cast<int>(node->rwgIndices.size());

			for (int r = 0; r < rwgNum; ++r) {
				const int rwgid = node->rwgIndices[r]->rwgid;
				const auto& V2 = Vfmj2[rwgid];

				const auto L2J = integrateL(
					V2, Bm_J2[nodeIdx], WGL_k2[0], WGL_phi_k2[0], tN2, pN2);
				const auto K2J = integrateK(
					V2, Bm_J2[nodeIdx], kp_lvl_k2[0], WGL_k2[0], WGL_phi_k2[0], tN2, pN2);
				const auto K2M = integrateK(
					V2, Bm_M2[nodeIdx], kp_lvl_k2[0], WGL_k2[0], WGL_phi_k2[0], tN2, pN2);
				const auto L2M = integrateL(
					V2, Bm_M2[nodeIdx], WGL_k2[0], WGL_phi_k2[0], tN2, pN2);

				Z_full[rwgNum][rwgid] += const2 * wave.eta2() * L2J;
				Z_full[rwgNum][N + rwgid] += -const2 * wave.eta1() * K2J;
				Z_full[rwgNum][rwgid] += const2 * wave.eta1() * K2M;
				Z_full[rwgNum][N + rwgid] += const2 * wave.eta2() * wave.epsilonR() * L2M;
			}
		}
	}

#pragma omp parallel for schedule(dynamic)
	for (int nodeIdx = 0; nodeIdx < nodenum; ++nodeIdx) {
		OCTree::Node* nd = octreeNodes_[maxLevel_][nodeIdx];
		for (int r = 0; r < static_cast<int>(nd->rwgIndices.size()); ++r) {
			int rwgid = nd->rwgIndices[r]->rwgid;

			for (int j = 0; j < static_cast<int>(Z_near[rwgid].size()); ++j) {
				Z_full[rwgid][Z_near_id[rwgid][j]] += x[Z_near_id[rwgid][j]] * Z_near[rwgid][j];
			}

			for (int j = 0; j < static_cast<int>(Z_near[N + rwgid].size()); ++j) {
				Z_full[N + rwgid][Z_near_id[N + rwgid][j]] += x[Z_near_id[N + rwgid][j]] * Z_near[N + rwgid][j];
			}
		}
	}
}