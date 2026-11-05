// MoM.cpp
#include "Zmartix.h"
#include "Vm.h"
#include "GMRES.h"
#include "CGS.h"
#include "RCS.h"
#include "MoM.h"

MoM::MoM(
	const RCSExportConfig& cfg, int selectMatrixSolver,
	const std::vector<RWGBase>& rwgs,
	const gaussPoints& gausspoint, const double E0, const EMSource& wave)
	: rwgs(rwgs), gausspoint(gausspoint), E0(E0), wave(wave), integralEquType_(static_cast<int>(cfg.integralEqu)), matrixSolverType_(selectMatrixSolver)
{
	row = static_cast<int>(rwgs.size());

	if (integralEquType_ == 0) {
		Zmartix::computeZ_pec_efie(rwgs, row, wave, gausspoint, Z_mom);

		//exportZ(cfg);

		if (cfg.monoDual == MonoDual::Dual) {
			mom_Dual_Pec_Efie(cfg);
		}

		if (cfg.monoDual == MonoDual::Mono) {
			mom_Mono_Pec_Efie(cfg);
		}
	}

	if (integralEquType_ == 1) {
		Zmartix::computeZ_pec_cfie(rwgs, row, wave, gausspoint, alpha, Z_mom);

		//exportZ(cfg);

		if (cfg.monoDual == MonoDual::Dual) {
			mom_Dual_Pec_Cfie(cfg);
		}

		if (cfg.monoDual == MonoDual::Mono) {
			mom_Mono_Pec_Cfie(cfg);
		}
	}

	if (integralEquType_ == 2) {
		Zmartix::computeZ_die_pmchwt(rwgs, row, wave, gausspoint, Z_mom);

		//exportZ(cfg);

		if (cfg.monoDual == MonoDual::Dual) {
			mom_Dual_Die_Pmchwt(cfg);
		}

		if (cfg.monoDual == MonoDual::Mono) {
			mom_Mono_Die_Pmchwt(cfg);
		}
	}
}

void MoM::exportZ(const RCSExportConfig& cfg) {
	auto zFile = cfg.openFile("_Zmatrix", "txt");

	zFile << std::setprecision(17);
	zFile << "row\tcol\treal\timag\n";

	for (int i = 0; i < static_cast<int>(Z_mom.size()); ++i) {
		for (int j = 0; j < static_cast<int>(Z_mom[i].size()); ++j) {
			zFile << i << '\t'
				<< j << '\t'
				<< Z_mom[i][j].real() << '\t'
				<< Z_mom[i][j].imag() << '\n';
		}
	}

	zFile.close();
	std::cout << "Z matrix export completed." << std::endl;
}

size_t MoM::computeMem() {
	size_t memZ_mom = my_memory::memory2D<std::complex<double>>(Z_mom);
	size_t memV = my_memory::memory1D<std::complex<double>>(Vm);
	return memZ_mom + memV;
}

void MoM::matrix_solver(
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
	bool converged = false;
	switch (matrixSolverType_)
	{
	case 0:
		std::cout << "Matrix solver: GMRES\n";
		converged = MGMRES::mgmres(n, x, rhs, itr_max, mr, tol_abs, tol_rel, ax);
		if (!converged) {
			std::cerr << "GMRES did not converge.\n";
		}
		break;
	case 1:
		std::cout << "Matrix solver: CGS\n";
		converged = MCGS::CGS::solve(n, x, rhs, itr_max, tol_rel, ax);
		if (!converged) {
			std::cerr << "CGS did not converge.\n";
		}
		break;
	default:
		throw std::invalid_argument(
			"Invalid matrixSolverType_ (must be 0 for GMRES or 1 for CGS)");
	}
}

void MoM::mom_Dual_Pec_Efie(const RCSExportConfig& cfg)
{
	// computeV lambda (EFIE)
	auto computeV = [](MoM& solver, double kInc[3], double eInc[3], double hInc[3]) {
		(void)hInc;
		RHS::computeV_E(solver.rwgs, solver.gausspoint,
			solver.wave.k1(), kInc, eInc, solver.Vm);
		};
	RCSUtils::computeDualStatic(*this, cfg, computeV);
}

void MoM::mom_Mono_Pec_Efie(const RCSExportConfig& cfg)
{
	auto computeV = [](MoM& solver, double kInc[3], double eInc[3], double hInc[3]) {
		(void)hInc;
		RHS::computeV_E(solver.rwgs, solver.gausspoint,
			solver.wave.k1(), kInc, eInc, solver.Vm);
		};
	RCSUtils::computeMonoStatic(*this, cfg, computeV);
}

void MoM::mom_Dual_Pec_Cfie(const RCSExportConfig& cfg)
{
	auto computeV = [](MoM& solver, double kInc[3], double eInc[3], double hInc[3]) {
		RHS::computeV_C(solver.rwgs, solver.gausspoint,
			solver.wave.k1(), kInc, eInc, hInc, alpha, solver.Vm);
		};
	RCSUtils::computeDualStatic(*this, cfg, computeV);
}

void MoM::mom_Mono_Pec_Cfie(const RCSExportConfig& cfg)
{
	auto computeV = [](MoM& solver, double kInc[3], double eInc[3], double hInc[3]) {
		RHS::computeV_C(solver.rwgs, solver.gausspoint,
			solver.wave.k1(), kInc, eInc, hInc, alpha, solver.Vm);
		};
	RCSUtils::computeMonoStatic(*this, cfg, computeV);
}

void MoM::mom_Dual_Die_Pmchwt(const RCSExportConfig& cfg)
{
	auto computeV = [](MoM& solver, double kInc[3], double eInc[3], double hInc[3]) {
		RHS::computeV_PMCHWT(solver.rwgs, solver.gausspoint,
			solver.wave.k1(), kInc, eInc, hInc, solver.Vm);
		};
	RCSUtils::computeDualStatic_PMCHWT(*this, cfg, computeV);
}

void MoM::mom_Mono_Die_Pmchwt(const RCSExportConfig& cfg)
{
	auto computeV = [](MoM& solver, double kInc[3], double eInc[3], double hInc[3]) {
		RHS::computeV_PMCHWT(solver.rwgs, solver.gausspoint,
			solver.wave.k1(), kInc, eInc, hInc, solver.Vm);
		};
	RCSUtils::computeMonoStatic_PMCHWT(*this, cfg, computeV);
}