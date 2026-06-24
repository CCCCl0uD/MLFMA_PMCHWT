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

	// k1: TF1
	AlgoMLFMM::computeTransferFactorMatrix(TF, octreeNodesDRvec_, kp_lvl_k1, maxLevel_, L_k1, wave.k1());

	// k1: Interpolation matrix
	AlgoMLFMM::computeInterpolationMatrix(interpol_k1, kp_lvl_k1, theta_level_k1, phi_level_k1, maxLevel_);

	initMLFMMStorage();

	if (integralEquType_ == 0) {
		// k1: Vsmi/Vfmj
		ProcessFMM::computeVsmi_Vfmj_Efie(Vsmi, Vfmj, octreeNodes_[maxLevel_],
			kp_lvl_k1[0], wave, row, gausspoint, alpha, 1);// Compute Vsmi and Vfmj (EFIE)

		// k1: Z_near (2N × 2N)
		Zmartix::computeZ_fmm_pec_efie(octreeNodes_[maxLevel_], rwgs, Z_near, Z_near_id,
			row, wave, gausspoint, 1);

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

		if (selectMono_Dual == "dual") {
			mlfmm_Dual_Pec_Cfie(cfg, pol_wave);
		}

		if (selectMono_Dual == "mono") {
			mlfmm_Mono_Pec_Cfie(cfg, pol_wave);
		}
	}

	if (integralEquType_ == 2) {
		// k1: Vsmi/Vfmj
		ProcessFMM::computeVsmi_Vfmj_Die(Vsmi, Vfmj,
			octreeNodes_[maxLevel_], kp_lvl_k1[0], wave.k1(), wave.eta1(), row, gausspoint,
			alpha, 1);

		// k2: Vsmi2/Vfmj2
		ProcessFMM::computeVsmi_Vfmj_Die(Vsmi2, Vfmj2,
			octreeNodes_[maxLevel_], kp_lvl_k2[0], wave.k2(), wave.eta2(), row, gausspoint,
			alpha, 1);

		// k2: TF2
		AlgoMLFMM::computeTransferFactorMatrix(TF2, octreeNodesDRvec_, kp_lvl_k2,
			maxLevel_, L_k2, wave.k2());

		// PMCHWT: Z_near (2N × 2N CSR)
		Zmartix::computeZ_fmm_die_pmchwt(octreeNodes_[maxLevel_], rwgs, Z_near, Z_near_id,
			row, wave, gausspoint, 1);

		if (selectMono_Dual == "dual") {
			mlfmm_Dual_Die_Pmchwt(cfg, pol_wave);
		}
		if (selectMono_Dual == "mono") {
			mlfmm_Mono_Die_Pmchwt(cfg, pol_wave);
		}
	}
}

void MLFMM::initMLFMMStorage() {
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
	// ---- k2 ----
	if (integralEquType_ == 2) {
		Sm2.resize(lvl);
		Bm2.resize(lvl);
		for (int level = 0; level < lvl; ++level) {
			int nodes = static_cast<int>(octreeNodes_[level + 2].size());
			int rev = levelSpan - 1 - level;
			int kp = static_cast<int>(kp_lvl_k2[rev].size());
			Sm2[level].resize(nodes);
			Bm2[level].resize(nodes);
			for (int node = 0; node < nodes; ++node) {
				Sm2[level][node].assign(kp, Complex3D());
				Bm2[level][node].assign(kp, Complex3D());
			}
		}
	}
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
		memSm_Bm += memory3D<Complex3D>(Sm2) + memory3D<Complex3D>(Bm2);
	}
	return memBase + memVsmi_Vfmi + memTF + memZ_near + memInterp + memSm_Bm + memV;
}

void MLFMM::matrix_solver(int n, std::complex<double> x[], std::complex<double> rhs[], int itr_max, int mr, double tol_abs, double tol_rel) {
	const int tN_k1 = static_cast<int>(theta_level_k1[0].size());
	const int pN_k1 = static_cast<int>(phi_level_k1[0].size());
	const int maxlvl = maxLevel_;

	auto ax = [this, tN_k1, pN_k1, maxlvl](int n, const std::complex<double>* x, std::complex<double>* w) {
		for (int i = 0; i < n; ++i) {
			w[i] = std::complex<double>(0.0, 0.0);
		}

		// ============ PEC ============
		if (integralEquType_ != 2) {
			ProcessMLFMM::clearSmBm(Sm, Bm);
			ProcessMLFMM::mexp(Sm, octreeNodes_[maxLevel_], Vsmi, kp_lvl_k1[0], x, 0, maxlvl);
			ProcessMLFMM::m2m(Sm, octreeNodes_, kp_lvl_k1, interpol_k1, maxlvl, levelSpan, wave.k1());
			ProcessMLFMM::m2l(Sm, Bm, octreeNodes_, TF, kp_lvl_k1, WGL_k1, WGL_phi_k1, phi_level_k1, maxlvl, levelSpan);
			ProcessMLFMM::l2l(Bm, octreeNodes_, kp_lvl_k1, phi_level_k1, interpol_k1, maxlvl, levelSpan, wave.k1());
			ProcessMLFMM::lexpPEC(w, octreeNodes_[maxLevel_], Vfmj, Bm, tN_k1, pN_k1, maxlvl,
				wave.k1(), wave.eta1(), Z_near, Z_near_id, x);
			return;
		}

		// ============ PMCHW ============
		const int N = row;
		const int tN_k2 = static_cast<int>(theta_level_k2[0].size());
		const int pN_k2 = static_cast<int>(phi_level_k2[0].size());
		const std::complex<double> const1 = wave.k1() * wave.k1() * wave.eta1() / (16.0 * Pi * Pi);
		const std::complex<double> const2 = wave.k2() * wave.k2() / (16.0 * Pi * Pi);

		// ---- Pass 1: k1 空间 J → L1(J)、K1(J) ----
		ProcessMLFMM::clearSmBm(Sm, Bm);
		ProcessMLFMM::mexp(Sm, octreeNodes_[maxLevel_], Vsmi, kp_lvl_k1[0], x, 0, maxlvl);
		ProcessMLFMM::m2m(Sm, octreeNodes_, kp_lvl_k1, interpol_k1, maxlvl, levelSpan, wave.k1());
		ProcessMLFMM::m2l(Sm, Bm, octreeNodes_, TF, kp_lvl_k1, WGL_k1, WGL_phi_k1, phi_level_k1, maxlvl, levelSpan);
		ProcessMLFMM::l2l(Bm, octreeNodes_, kp_lvl_k1, phi_level_k1, interpol_k1, maxlvl, levelSpan, wave.k1());
		ProcessMLFMM::lexpPMCHWT(w, const1, -const1, 0, N, N,
			octreeNodes_[maxLevel_], Vfmj, Bm, kp_lvl_k1[0], tN_k1, pN_k1, maxlvl);

		// ---- Pass 2: k1 空间 M → K1(M)、L1(M) ----
		ProcessMLFMM::clearSmBm(Sm, Bm);
		ProcessMLFMM::mexp(Sm, octreeNodes_[maxLevel_], Vsmi, kp_lvl_k1[0], x, N, maxlvl);
		ProcessMLFMM::m2m(Sm, octreeNodes_, kp_lvl_k1, interpol_k1, maxlvl, levelSpan, wave.k1());
		ProcessMLFMM::m2l(Sm, Bm, octreeNodes_, TF, kp_lvl_k1, WGL_k1, WGL_phi_k1, phi_level_k1, maxlvl, levelSpan);
		ProcessMLFMM::l2l(Bm, octreeNodes_, kp_lvl_k1, phi_level_k1, interpol_k1, maxlvl, levelSpan, wave.k1());
		ProcessMLFMM::lexpPMCHWT(w, const1, const1, N, 0, N,
			octreeNodes_[maxLevel_], Vfmj, Bm, kp_lvl_k1[0], tN_k1, pN_k1, maxlvl);

		// ---- Pass 3: k2 空间 J → L2(J)、K2(J) ----
		ProcessMLFMM::clearSmBm(Sm2, Bm2);
		ProcessMLFMM::mexp(Sm2, octreeNodes_[maxLevel_], Vsmi2, kp_lvl_k2[0], x, 0, maxlvl);
		ProcessMLFMM::m2m(Sm2, octreeNodes_, kp_lvl_k2, interpol_k2, maxlvl, levelSpan, wave.k2());
		ProcessMLFMM::m2l(Sm2, Bm2, octreeNodes_, TF2, kp_lvl_k2, WGL_k2, WGL_phi_k2, phi_level_k2, maxlvl, levelSpan);
		ProcessMLFMM::l2l(Bm2, octreeNodes_, kp_lvl_k2, phi_level_k2, interpol_k2, maxlvl, levelSpan, wave.k2());
		ProcessMLFMM::lexpPMCHWT(w, const2 * wave.eta2(), -const2 * wave.eta1(), 0, N, N,
			octreeNodes_[maxLevel_], Vfmj2, Bm2, kp_lvl_k2[0], tN_k2, pN_k2, maxlvl);

		// ---- Pass 4: k2 空间 M → K2(M)、L2(M) ----
		ProcessMLFMM::clearSmBm(Sm2, Bm2);
		ProcessMLFMM::mexp(Sm2, octreeNodes_[maxLevel_], Vsmi2, kp_lvl_k2[0], x, N, maxlvl);
		ProcessMLFMM::m2m(Sm2, octreeNodes_, kp_lvl_k2, interpol_k2, maxlvl, levelSpan, wave.k2());
		ProcessMLFMM::m2l(Sm2, Bm2, octreeNodes_, TF2, kp_lvl_k2, WGL_k2, WGL_phi_k2, phi_level_k2, maxlvl, levelSpan);
		ProcessMLFMM::l2l(Bm2, octreeNodes_, kp_lvl_k2, phi_level_k2, interpol_k2, maxlvl, levelSpan, wave.k2());
		ProcessMLFMM::lexpPMCHWT(w, const2 * wave.eta2() * wave.epsilonR(), const2 * wave.eta1(), N, 0, N,
			octreeNodes_[maxLevel_], Vfmj2, Bm2, kp_lvl_k2[0], tN_k2, pN_k2, maxlvl);

		// ---- 近场: 2N×2N ----
		const int finestN = static_cast<int>(octreeNodes_[maxLevel_].size());
#pragma omp parallel for schedule(dynamic)
		for (int nodeIdx = 0; nodeIdx < finestN; ++nodeIdx) {
			OCTree::Node* nd = octreeNodes_[maxLevel_][nodeIdx];
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
		};
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
			solver.wave.k1(), solver.wave.k2(), kInc, eInc, hInc, solver.Vm);
		};
	RCSUtils::computeDualStatic_PMCHWT(*this, cfg, computeV);
}

void MLFMM::mlfmm_Mono_Die_Pmchwt(const RCSExportConfig& cfg, const std::string pol_wave)
{
	auto computeV = [](MLFMM& solver, double kInc[3], double eInc[3], double hInc[3]) {
		RHS::computeV_PMCHWT(solver.rwgs, solver.gausspoint,
			solver.wave.k1(), solver.wave.k2(), kInc, eInc, hInc, solver.Vm);
		};
	RCSUtils::computeMonoStatic_PMCHWT(*this, cfg, computeV);
}