// MLFMM.cpp
#include "MLFMM.h"
#include "AlgoFMM.h"
#include "AlgoMLFMM.h"
#include "GMRES.h"
#include "CGS.h"
#include "ProcessFMM.h"
#include "ProcessMLFMM.h"
#include "Zmartix.h"
#include "Vm.h"
#include "RCS.h"

/**************************************************************************************************/
/**************************************************************************************************/
MLFMM::MLFMM(
	const RCSExportConfig& cfg, const int selectIntegralEqu, const int selectMatrixSolver, const std::string selectMono_Dual, const std::string pol_wave,
	const std::vector<std::vector<OCTree::Node*>>& octreeNodes,
	const std::vector<std::vector<OCTree::nodePoint>>& octreeNodesDRvec_,
	const std::vector<RWGBase>& rwgs, const int maxLevel_,
	const gaussPoints& gausspoint, const double E0, const EMSource& wave)
	: octreeNodes_(octreeNodes), octreeNodesDRvec_(octreeNodesDRvec_), rwgs(rwgs), maxLevel_(maxLevel_), gausspoint(gausspoint), E0(E0), wave(wave), integralEquType_(selectIntegralEqu), matrixSolverType_(selectMatrixSolver)
{
	// Compute base data: k, top level theta, L, levelSpan
	// Compute Vsmi/Vfmj matrices theta/phi angles and weights from maxLevel_ down to level 2
	AlgoFMM::computeBase(octreeNodes_, wave, integralEquType_, maxLevel_, static_cast<int>(rwgs.size()), L_k1, L_k2, row, levelSpan,
		WGL_k1, WGL_k2, WGL_phi_k1, WGL_phi_k2, theta_level_k1, theta_level_k2, phi_level_k1, phi_level_k2, kp_lvl_k1, kp_lvl_k2);

	if (integralEquType_ != 2) {
		// k1: TF1
		AlgoMLFMM::computeTransferFactorMatrix(TF, octreeNodesDRvec_, kp_lvl_k1, maxLevel_, L_k1, wave.k1());

		// k1: Interpolation matrix
		AlgoMLFMM::computeInterpolationMatrix(interpol_k1, kp_lvl_k1, theta_level_k1, phi_level_k1, maxLevel_);

		if (integralEquType_ == 0) {
			// k1: Vsmi/Vfmj
			ProcessFMM::computeVsmi_Vfmj_Efie(Vsmi, Vfmj, octreeNodes_[maxLevel_],
				kp_lvl_k1[0], wave, row, gausspoint, alpha, 1);// Compute Vsmi and Vfmj (EFIE)

			// k1: Z_near (2N × 2N)
			Zmartix::computeZ_fmm_pec_efie(octreeNodes_[maxLevel_], rwgs, Z_near, Z_near_id,
				row, wave, gausspoint, 1);

			initPECWorkspace();

			if (selectMono_Dual == "dual") {
				mlfmm_Dual_Pec_Efie(cfg, pol_wave);
			}
			if (selectMono_Dual == "mono") {
				mlfmm_Mono_Pec_Efie(cfg, pol_wave);
			}
		}

		if (integralEquType_ == 1) {
			// k1: Vsmi/Vfmj
			ProcessFMM::computeVsmi_Vfmj_Cfie(Vsmi, Vfmj, octreeNodes_[maxLevel_],
				kp_lvl_k1[0], wave, row, gausspoint, alpha, 1);// Compute Vsmi and Vfmj (CFIE)

			// k1: Z_near (2N × 2N)
			Zmartix::computeZ_fmm_pec_cfie(octreeNodes_[maxLevel_], rwgs, Z_near, Z_near_id,
				row, wave, gausspoint, alpha, 1);

			initPECWorkspace();

			if (selectMono_Dual == "dual") {
				mlfmm_Dual_Pec_Cfie(cfg, pol_wave);
			}
			if (selectMono_Dual == "mono") {
				mlfmm_Mono_Pec_Cfie(cfg, pol_wave);
			}
		}
	}
	else {
		// k1: Vsmi/Vfmj
		ProcessFMM::computeVsmi_Vfmj_Die(Vsmi, Vfmj,
			octreeNodes_[maxLevel_], kp_lvl_k1[0], wave.k1(), row, gausspoint,
			alpha, 1);

		// k2: Vsmi2/Vfmj2
		ProcessFMM::computeVsmi_Vfmj_Die(Vsmi2, Vfmj2,
			octreeNodes_[maxLevel_], kp_lvl_k2[0], wave.k2(), row, gausspoint,
			alpha, 1);

		// k1: TF1
		AlgoMLFMM::computeTransferFactorMatrix(TF, octreeNodesDRvec_, kp_lvl_k1,
			maxLevel_, L_k1, wave.k1());

		// k2: TF2
		AlgoMLFMM::computeTransferFactorMatrix(TF2, octreeNodesDRvec_, kp_lvl_k2,
			maxLevel_, L_k2, wave.k2());

		// k1: Interpolation matrix
		AlgoMLFMM::computeInterpolationMatrix(interpol_k1, kp_lvl_k1,
			theta_level_k1, phi_level_k1, maxLevel_);

		// k2: Interpolation matrix
		AlgoMLFMM::computeInterpolationMatrix(interpol_k2, kp_lvl_k2,
			theta_level_k2, phi_level_k2, maxLevel_);

		// PMCHWT: Z_near (2N × 2N CSR)
		Zmartix::computeZ_fmm_die_pmchwt(octreeNodes_[maxLevel_], rwgs, Z_near, Z_near_id,
			row, wave, gausspoint, 1);

		initDIEWorkspace();

		if (selectMono_Dual == "dual") {
			mlfmm_Dual_Die_Pmchwt(cfg, pol_wave);
		}
		if (selectMono_Dual == "mono") {
			mlfmm_Mono_Die_Pmchwt(cfg, pol_wave);
		}
	}
}

void MLFMM::initPECWorkspace() {
	int lvl = maxLevel_ - 1;
	// ---- k1 ----
	Sm.resize(lvl);
	Bm.resize(lvl);
	for (int level = 0; level < lvl; ++level) {
		int nodes = static_cast<int>(octreeNodes_[level + 2].size());
		int rev = levelSpan - 1 - level;
		int kp = static_cast<int>(kp_lvl_k1[rev].size());
		Sm[level].resize(nodes);
		Bm[level].resize(nodes);
		for (int node = 0; node < nodes; ++node) {
			Sm[level][node].assign(kp, Complex3D());
			Bm[level][node].assign(kp, Complex3D());
		}
	}
}

void MLFMM::initDIEWorkspace() {
	const int lvl = maxLevel_ - 1;

	auto initWork = [&](std::vector<std::vector<std::vector<Complex3D>>>& SmWork,
		std::vector<std::vector<std::vector<Complex3D>>>& BmWork,
		const std::vector<std::vector<kp_Point>>& kpLevel) {
			SmWork.resize(lvl);
			BmWork.resize(lvl);

			for (int level = 0; level < lvl; ++level) {
				const int nodes = static_cast<int>(octreeNodes_[level + 2].size());
				const int rev = levelSpan - 1 - level;
				const int kp = static_cast<int>(kpLevel[rev].size());

				SmWork[level].resize(nodes);
				BmWork[level].resize(nodes);

				for (int node = 0; node < nodes; ++node) {
					SmWork[level][node].assign(kp, Complex3D{});
					BmWork[level][node].assign(kp, Complex3D{});
				}
			}
		};

	initWork(Sm_J1, Bm_J1, kp_lvl_k1);
	initWork(Sm_M1, Bm_M1, kp_lvl_k1);
	initWork(Sm_J2, Bm_J2, kp_lvl_k2);
	initWork(Sm_M2, Bm_M2, kp_lvl_k2);
}

size_t MLFMM::computeMem() {
	size_t memBase = memory2D<kp_Point>(kp_lvl_k1) + memory2D<double>(theta_level_k1) + memory2D<double>(phi_level_k1) +
		memory2D<double>(WGL_k1) + memory1D<double>(WGL_phi_k1);
	size_t memVsmi_Vfmi = memory2D<Complex3D>(Vsmi) + memory2D<Complex3D>(Vfmj);
	size_t memTF = memory3D<std::complex<double>>(TF);
	size_t memZ_near = memory2D<std::complex<double>>(Z_near);
	size_t memInterp = memory2D<InterpData>(interpol_k1);
	size_t memSm_Bm = memory3D<Complex3D>(Sm) + memory3D<Complex3D>(Bm);
	size_t memV = memory1D<std::complex<double>>(Vm);
	if (integralEquType_ == 2) {
		memBase += memory2D<kp_Point>(kp_lvl_k2) + memory2D<double>(theta_level_k2) + memory2D<double>(phi_level_k2) +
			memory2D<double>(WGL_k2) + memory1D<double>(WGL_phi_k2);
		memVsmi_Vfmi += memory2D<Complex3D>(Vsmi2) + memory2D<Complex3D>(Vfmj2);
		memTF += memory3D<std::complex<double>>(TF2);
		memInterp += memory2D<InterpData>(interpol_k2);
		memSm_Bm += memory3D<Complex3D>(Sm_J1) + memory3D<Complex3D>(Bm_J1)
			+ memory3D<Complex3D>(Sm_M1) + memory3D<Complex3D>(Bm_M1)
			+ memory3D<Complex3D>(Sm_J2) + memory3D<Complex3D>(Bm_J2)
			+ memory3D<Complex3D>(Sm_M2) + memory3D<Complex3D>(Bm_M2);
	}
	return memBase + memVsmi_Vfmi + memTF + memZ_near + memInterp + memSm_Bm + memV;
}

void MLFMM::matrix_solver(int n, std::complex<double> x[], std::complex<double> rhs[], int itr_max, int mr, double tol_abs, double tol_rel) {
	const int EquType_ = integralEquType_;

	auto ax = [this, EquType_](int n, const std::complex<double>* x, std::complex<double>* w) {
		for (int i = 0; i < n; ++i) {
			w[i] = std::complex<double>(0.0, 0.0);
		}

		// ============ PEC ============
		if (EquType_ != 2) {
			const int maxlvl = maxLevel_;

			const int tN_k1 = static_cast<int>(theta_level_k1[0].size());
			const int pN_k1 = static_cast<int>(phi_level_k1[0].size());
			const int nodeNum = static_cast<int>(octreeNodes_[maxlvl].size());

			const std::complex<double> k_ = wave.k1();
			const std::complex<double> eta_ = wave.eta1();

			ProcessMLFMM::clearSmBm(Sm, Bm);
			ProcessMLFMM::mexp(Sm, octreeNodes_[maxlvl], Vsmi, kp_lvl_k1[0], x, 0, maxlvl);
			ProcessMLFMM::m2m(Sm, octreeNodes_, kp_lvl_k1, interpol_k1, maxlvl, levelSpan, wave.k1());
			ProcessMLFMM::m2l(Sm, Bm, octreeNodes_, TF, kp_lvl_k1, WGL_k1, WGL_phi_k1, phi_level_k1, maxlvl, levelSpan);
			ProcessMLFMM::l2l(Bm, octreeNodes_, kp_lvl_k1, phi_level_k1, interpol_k1, maxlvl, levelSpan, wave.k1());

#pragma omp parallel for schedule(dynamic)
			for (int nodeIdx = 0; nodeIdx < nodeNum; ++nodeIdx) {
				OCTree::Node* node = octreeNodes_[maxlvl][nodeIdx];
				const int rwgNum = static_cast<int>(node->rwgIndices.size());
				const auto& B = Bm[maxlvl - 2][nodeIdx];

				for (int r = 0; r < rwgNum; ++r) {
					int rwgid = node->rwgIndices[r]->rwgid;
					const auto& V = Vfmj[rwgid];
					std::complex<double> val(0.0, 0.0);
					for (int tt = 0; tt < tN_k1; ++tt) {
						for (int pp = 0; pp < pN_k1; ++pp) {
							val += ProcessMLFMM::lexpKpDot(B[tt * pN_k1 + pp], V[tt * pN_k1 + pp]);
						}
					}
					val *= (k_ * eta_ / (4.0 * Pi)) * (k_ / (4.0 * Pi));

					const int near = static_cast<int>(Z_near[rwgid].size());
					for (int j = 0; j < near; ++j) {
						val += x[Z_near_id[rwgid][j]] * Z_near[rwgid][j];
					}
					w[rwgid] = val;
				}
			}
			return;
		}
		else {
			// ===== PMCHWT matrix-vector product =====
			// Unknown ordering:
			//   x[0:N)     = \eta_0 * J coefficients
			//   x[N:2*N)   = M coefficients
			//
			// Equation ordering:
			//   w[0:N)     = electric-current testing block
			//   w[N:2*N)   = magnetic-current testing block
			const int N = row;
			const int maxlvl = maxLevel_;
			const int finestIdx = maxlvl - 2;

			const std::complex<double> etai = wave.etai();

			const std::complex<double> const1 = wave.k1() * wave.k1() / (16.0 * Pi * Pi);
			const std::complex<double> const2 = wave.k2() * wave.k2() / (16.0 * Pi * Pi);

			ProcessMLFMM::clearSmBm(Sm_J1, Bm_J1);
			ProcessMLFMM::clearSmBm(Sm_M1, Bm_M1);
			ProcessMLFMM::clearSmBm(Sm_J2, Bm_J2);
			ProcessMLFMM::clearSmBm(Sm_M2, Bm_M2);

			// ---- Pass 1: k1, J → L1(J)、K1(J) ----
			ProcessMLFMM::mexp(Sm_J1, octreeNodes_[maxlvl], Vsmi, kp_lvl_k1[0], x, 0, maxlvl);
			ProcessMLFMM::m2m(Sm_J1, octreeNodes_, kp_lvl_k1, interpol_k1, maxlvl, levelSpan, wave.k1());
			ProcessMLFMM::m2l(Sm_J1, Bm_J1, octreeNodes_, TF, kp_lvl_k1, WGL_k1, WGL_phi_k1, phi_level_k1, maxlvl, levelSpan);
			ProcessMLFMM::l2l(Bm_J1, octreeNodes_, kp_lvl_k1, phi_level_k1, interpol_k1, maxlvl, levelSpan, wave.k1());

			// ---- Pass 2: k1, M → K1(M)、L1(M) ----
			ProcessMLFMM::mexp(Sm_M1, octreeNodes_[maxlvl], Vsmi, kp_lvl_k1[0], x, N, maxlvl);
			ProcessMLFMM::m2m(Sm_M1, octreeNodes_, kp_lvl_k1, interpol_k1, maxlvl, levelSpan, wave.k1());
			ProcessMLFMM::m2l(Sm_M1, Bm_M1, octreeNodes_, TF, kp_lvl_k1, WGL_k1, WGL_phi_k1, phi_level_k1, maxlvl, levelSpan);
			ProcessMLFMM::l2l(Bm_M1, octreeNodes_, kp_lvl_k1, phi_level_k1, interpol_k1, maxlvl, levelSpan, wave.k1());

			// ---- Pass 3: k2, J → L2(J)、K2(J) ----
			ProcessMLFMM::mexp(Sm_J2, octreeNodes_[maxlvl], Vsmi2, kp_lvl_k2[0], x, 0, maxlvl);
			ProcessMLFMM::m2m(Sm_J2, octreeNodes_, kp_lvl_k2, interpol_k2, maxlvl, levelSpan, wave.k2());
			ProcessMLFMM::m2l(Sm_J2, Bm_J2, octreeNodes_, TF2, kp_lvl_k2, WGL_k2, WGL_phi_k2, phi_level_k2, maxlvl, levelSpan);
			ProcessMLFMM::l2l(Bm_J2, octreeNodes_, kp_lvl_k2, phi_level_k2, interpol_k2, maxlvl, levelSpan, wave.k2());

			// ---- Pass 4: k2, M → K2(M)、L2(M) ----
			ProcessMLFMM::mexp(Sm_M2, octreeNodes_[maxlvl], Vsmi2, kp_lvl_k2[0], x, N, maxlvl);
			ProcessMLFMM::m2m(Sm_M2, octreeNodes_, kp_lvl_k2, interpol_k2, maxlvl, levelSpan, wave.k2());
			ProcessMLFMM::m2l(Sm_M2, Bm_M2, octreeNodes_, TF2, kp_lvl_k2, WGL_k2, WGL_phi_k2, phi_level_k2, maxlvl, levelSpan);
			ProcessMLFMM::l2l(Bm_M2, octreeNodes_, kp_lvl_k2, phi_level_k2, interpol_k2, maxlvl, levelSpan, wave.k2());

			//---- far field ----
			auto integrateL = [](const std::vector<Complex3D>& V,
				const std::vector<Complex3D>& B) -> std::complex<double> {
					std::complex<double> val(0.0, 0.0);
					const int K = static_cast<int>(V.size());
					for (int idx = 0; idx < K; ++idx) {
						val += dot(V[idx], B[idx]);
					}
					return val;
				};

			auto integrateK = [](const std::vector<Complex3D>& V,
				const std::vector<Complex3D>& B,
				const std::vector<kp_Point>& kp) -> std::complex<double> {
					std::complex<double> val(0.0, 0.0);
					const int K = static_cast<int>(V.size());
					for (int idx = 0; idx < K; ++idx) {
						const double khat[3] = { -kp[idx].x, -kp[idx].y, -kp[idx].z };

						Complex3D kCrossV;
						cross(kCrossV, khat, V[idx]);

						val += dot(kCrossV, B[idx]);
					}
					return val;
				};

			// ---- far field: finest-level disaggregation ----
			const int finestN = static_cast<int>(octreeNodes_[maxlvl].size());

#pragma omp parallel for schedule(dynamic)
			for (int nodeIdx = 0; nodeIdx < finestN; ++nodeIdx) {
				OCTree::Node* nd = octreeNodes_[maxlvl][nodeIdx];

				const auto& BJ1 = Bm_J1[finestIdx][nodeIdx];
				const auto& BM1 = Bm_M1[finestIdx][nodeIdx];
				const auto& BJ2 = Bm_J2[finestIdx][nodeIdx];
				const auto& BM2 = Bm_M2[finestIdx][nodeIdx];

				for (int ri = 0; ri < static_cast<int>(nd->rwgIndices.size()); ++ri) {
					const int rwgid = nd->rwgIndices[ri]->rwgid;

					const auto& V1 = Vfmj[rwgid];
					const auto& V2 = Vfmj2[rwgid];

					const auto L1J = integrateL(V1, BJ1);
					const auto K1J = integrateK(V1, BJ1, kp_lvl_k1[0]);
					const auto K1M = integrateK(V1, BM1, kp_lvl_k1[0]);
					const auto L1M = integrateL(V1, BM1);

					const auto L2J = integrateL(V2, BJ2);
					const auto K2J = integrateK(V2, BJ2, kp_lvl_k2[0]);
					const auto K2M = integrateK(V2, BM2, kp_lvl_k2[0]);
					const auto L2M = integrateL(V2, BM2);

					w[rwgid] += const1 * L1J;
					w[rwgid] += -const1 * K1M;

					w[N + rwgid] += const1 * K1J;
					w[N + rwgid] += const1 * L1M;

					w[rwgid] += const2 * etai * L2J;
					w[rwgid] += -const2 * K2M;

					w[N + rwgid] += const2 * K2J;
					w[N + rwgid] += const2 * L2M / etai;
				}
			}

			// ---- near field ----
#pragma omp parallel for schedule(dynamic)
			for (int nodeIdx = 0; nodeIdx < finestN; ++nodeIdx) {
				OCTree::Node* nd = octreeNodes_[maxlvl][nodeIdx];
				for (int ri = 0; ri < static_cast<int>(nd->rwgIndices.size()); ++ri) {
					int rwgid = nd->rwgIndices[ri]->rwgid;
					for (int j = 0; j < static_cast<int>(Z_near[rwgid].size()); ++j) {
						w[rwgid] += x[Z_near_id[rwgid][j]] * Z_near[rwgid][j];
					}
					for (int j = 0; j < static_cast<int>(Z_near[N + rwgid].size()); ++j) {
						w[N + rwgid] += x[Z_near_id[N + rwgid][j]] * Z_near[N + rwgid][j];
					}
				}
			}
			return;
		}
		};
	bool converged = false;

	switch (matrixSolverType_) {
	case 0:
		std::cout << "Matrix solver: GMRES (MLFMM)\n";
		converged = MGMRES::mgmres(
			n, x, rhs, itr_max, mr, tol_abs, tol_rel, ax);
		if (!converged) {
			std::cerr << "GMRES did not converge (MLFMM).\n";
		}
		break;

	case 1:
		std::cout << "Matrix solver: CGS (MLFMM)\n";
		converged = MCGS::CGS::solve(
			n, x, rhs, itr_max, tol_rel, ax);
		if (!converged) {
			std::cerr
				<< "CGS did not converge or encountered breakdown (MLFMM).\n";
		}
		break;

	default:
		throw std::invalid_argument(
			"Invalid matrixSolverType_ in MLFMM "
			"(must be 0 for GMRES or 1 for CGS)");
	}
}

void MLFMM::mlfmm_Dual_Pec_Efie(const RCSExportConfig& cfg, const std::string pol_wave)
{
	auto computeV = [](MLFMM& solver, double kInc[3], double eInc[3], double hInc[3]) {
		(void)hInc;
		RHS::computeV_E(solver.rwgs, solver.gausspoint,
			solver.wave.k1(), kInc, eInc, solver.Vm);
		};
	RCSUtils::computeDualStatic(*this, cfg, computeV);
}

void MLFMM::mlfmm_Mono_Pec_Efie(const RCSExportConfig& cfg, const std::string pol_wave)
{
	auto computeV = [](MLFMM& solver, double kInc[3], double eInc[3], double hInc[3]) {
		(void)hInc;
		RHS::computeV_E(solver.rwgs, solver.gausspoint,
			solver.wave.k1(), kInc, eInc, solver.Vm);
		};
	RCSUtils::computeMonoStatic(*this, cfg, computeV);
}

void MLFMM::mlfmm_Dual_Pec_Cfie(const RCSExportConfig& cfg, const std::string pol_wave)
{
	auto computeV = [](MLFMM& solver, double kInc[3], double eInc[3], double hInc[3]) {
		RHS::computeV_C(solver.rwgs, solver.gausspoint,
			solver.wave.k1(), kInc, eInc, hInc, alpha, solver.Vm);
		};
	RCSUtils::computeDualStatic(*this, cfg, computeV);
}

void MLFMM::mlfmm_Mono_Pec_Cfie(const RCSExportConfig& cfg, const std::string pol_wave)
{
	auto computeV = [this](MLFMM& solver, double kInc[3], double eInc[3], double hInc[3]) {
		RHS::computeV_C(solver.rwgs, solver.gausspoint,
			solver.wave.k1(), kInc, eInc, hInc, alpha, solver.Vm);
		};
	RCSUtils::computeMonoStatic(*this, cfg, computeV);
}

void MLFMM::mlfmm_Dual_Die_Pmchwt(const RCSExportConfig& cfg, const std::string pol_wave) {
	auto computeV = [](MLFMM& solver, double kInc[3], double eInc[3], double hInc[3]) {
		RHS::computeV_PMCHWT(solver.rwgs, solver.gausspoint,
			solver.wave.k1(), kInc, eInc, hInc, solver.Vm);
		};
	RCSUtils::computeDualStatic_PMCHWT(*this, cfg, computeV);
}

void MLFMM::mlfmm_Mono_Die_Pmchwt(const RCSExportConfig& cfg, const std::string pol_wave)
{
	auto computeV = [](MLFMM& solver, double kInc[3], double eInc[3], double hInc[3]) {
		RHS::computeV_PMCHWT(solver.rwgs, solver.gausspoint,
			solver.wave.k1(), kInc, eInc, hInc, solver.Vm);
		};
	RCSUtils::computeMonoStatic_PMCHWT(*this, cfg, computeV);
}