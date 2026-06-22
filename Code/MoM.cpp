// MoM.cpp
#include "Zmartix.h"
#include "Vm.h"
#include "GMRES.h"
#include "RCS.h"
#include "MoM.h"

MoM::MoM(
	const RCSExportConfig& cfg, const int selectIntegralEqu, const std::string selectMono_Dual, const std::string pol_wave,
	const std::vector<std::vector<OCTree::Node*>>& octreeNodes, const std::vector<RWGBase>& rwgs, const int maxLevel_,
	const gaussPoints& gausspoint, const double E0, const EMSource& wave)
	: octreeNodes_(octreeNodes), rwgs(rwgs), maxLevel_(maxLevel_), gausspoint(gausspoint), E0(E0), wave(wave), integralEquType_(selectIntegralEqu)
{
	row = static_cast<int>(rwgs.size());

	if (integralEquType_ == 0) {
		Zmartix::computeZ_pec_efie(rwgs, row, wave, gausspoint, Z_mom);

		//exportZ(cfg);

		if (selectMono_Dual == "dual") {
			mom_Dual_Pec_Efie(cfg, pol_wave);
		}

		if (selectMono_Dual == "mono") {
			mom_Mono_Pec_Efie(cfg, pol_wave);
		}
	}

	if (integralEquType_ == 1) {
		Zmartix::computeZ_pec_cfie(rwgs, row, wave, gausspoint, alpha, Z_mom);

		//exportZ(cfg);

		if (selectMono_Dual == "dual") {
			mom_Dual_Pec_Cfie(cfg, pol_wave);
		}

		if (selectMono_Dual == "mono") {
			mom_Mono_Pec_Cfie(cfg, pol_wave);
		}
	}

	if (integralEquType_ == 2) {
		Zmartix::computeZ_die_pmchwt(rwgs, row, wave, gausspoint, Z_mom);

		//exportZ(cfg);

		if (selectMono_Dual == "dual") {
			mom_Dual_Die_Pmchwt(cfg, pol_wave);
		}

		if (selectMono_Dual == "mono") {
			mom_Mono_Die_Pmchwt(cfg, pol_wave);
		}
	}
}

void MoM::exportZ(const RCSExportConfig& cfg) {
	auto zFile = cfg.openFile("_Zmatrix", "txt");
	if (!zFile.is_open()) {
		std::cerr << "Error: Failed to open file for Z matrix export." << std::endl;
		return;
	}

	zFile << "Z_mom Matrix\n";

	for (int i = 0; i < row; i++)
	{
		for (int j = 0; j < Z_mom[i].size(); j++)
		{
			zFile << i << "\t" << j << "\t" << Z_mom[i][j] << "\n";
		}

		zFile << "\n";
	}
	zFile.close();
	std::cout << "Z matrix export completed." << std::endl;
}

size_t MoM::computeMem() {
	size_t memZ_mom = memory2D<std::complex<double>>(Z_mom);
	size_t memV = memory1D<std::complex<double>>(Vm);
	return memZ_mom + memV;
}

void MoM::mgmres_solver(
	int n, std::complex<double> x[], std::complex<double> rhs[],
	int itr_max, int mr, double tol_abs, double tol_rel)
{
	auto ax = [this](int n, const std::complex<double>* x, std::complex<double>* w) {
#pragma omp parallel for
		for (int i = 0; i < n; ++i) {
			std::complex<double> sum = 0.0;
			for (int j = 0; j < n; ++j) {
				sum += Z_mom[i][j] * x[j];
			}
			w[i] = sum;
		}
		};
	bool converged = MGMRES::mgmres(n, x, rhs, itr_max, mr, tol_abs, tol_rel, ax);
	if (!converged) {
		std::cerr << "MGMRES did not converge.\n";
	}
}

void MoM::mom_Dual_Pec_Efie(const RCSExportConfig& cfg, const std::string pol_wave)
{
	// computeV lambda (EFIE)
	auto computeV = [](MoM& solver, double kInc[3], double eInc[3], double hInc[3]) {
		(void)hInc;
		RHS::computeV_E(solver.rwgs, solver.gausspoint,
			solver.wave.k1(), kInc, eInc, solver.Vm);
		};
	RCSUtils::computeDualStatic(*this, cfg, computeV);
}

void MoM::mom_Mono_Pec_Efie(const RCSExportConfig& cfg, const std::string pol_wave)
{
	auto computeV = [](MoM& solver, double kInc[3], double eInc[3], double hInc[3]) {
		(void)hInc;
		RHS::computeV_E(solver.rwgs, solver.gausspoint,
			solver.wave.k1(), kInc, eInc, solver.Vm);
		};
	RCSUtils::computeMonoStatic(*this, cfg, computeV);
}

void MoM::mom_Dual_Pec_Cfie(const RCSExportConfig& cfg, const std::string pol_wave)
{
	auto computeV = [](MoM& solver, double kInc[3], double eInc[3], double hInc[3]) {
		RHS::computeV_C(solver.rwgs, solver.gausspoint,
			solver.wave.k1(), kInc, eInc, hInc, alpha, solver.Vm);
		};
	RCSUtils::computeDualStatic(*this, cfg, computeV);
}

void MoM::mom_Mono_Pec_Cfie(const RCSExportConfig& cfg, const std::string pol_wave)
{
	auto computeV = [](MoM& solver, double kInc[3], double eInc[3], double hInc[3]) {
		RHS::computeV_C(solver.rwgs, solver.gausspoint,
			solver.wave.k1(), kInc, eInc, hInc, alpha, solver.Vm);
		};
	RCSUtils::computeMonoStatic(*this, cfg, computeV);
}

void MoM::mom_Dual_Die_Pmchwt(const RCSExportConfig& cfg, const std::string pol_wave)
{
	auto computeV = [](MoM& solver, double kInc[3], double eInc[3], double hInc[3]) {
		RHS::computeV_PMCHWT(solver.rwgs, solver.gausspoint,
			solver.wave.k1(), solver.wave.eta1(), kInc, eInc, hInc, solver.Vm);
		};
	RCSUtils::computeDualStatic_PMCHWT(*this, cfg, computeV);
}

void MoM::mom_Mono_Die_Pmchwt(const RCSExportConfig& cfg, const std::string pol_wave)
{
	auto computeV = [](MoM& solver, double kInc[3], double eInc[3], double hInc[3]) {
		RHS::computeV_PMCHWT(solver.rwgs, solver.gausspoint,
			solver.wave.k1(), solver.wave.eta1(), kInc, eInc, hInc, solver.Vm);
		};
	RCSUtils::computeMonoStatic_PMCHWT(*this, cfg, computeV);
}