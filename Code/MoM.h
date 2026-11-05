// MoM.h
#pragma once
#ifndef MOM_H
#define MOM_H

#include <vector>

#include "MyStruct.h"
#include "RCSExportConfig.h"
#include "EMSource.h"
#include "GaussPoints.h"

class MoM {
public:
	MoM(const RCSExportConfig& cfg, int selectMatrixSolver,
		const std::vector<RWGBase>& rwgs,
		const gaussPoints& gausspoint, const double E0, const EMSource& wave);

	const std::vector<RWGBase>& rwgs;
	const double E0;
	int row;                                              // RWG基函数个数
	const gaussPoints& gausspoint;
	const EMSource& wave;
	const int integralEquType_;
	const int matrixSolverType_;
	std::vector<std::complex<double>> Vm;                 // 右边向量
	std::vector<std::vector<std::complex<double>>> Z_mom;

	size_t computeMem();

	void matrix_solver(int n, std::complex<double> x[], std::complex<double> rhs[], int itr_max, int mr, double tol_abs, double tol_rel);

private:
	void mom_Dual_Pec_Efie(const RCSExportConfig& cfg);
	void mom_Mono_Pec_Efie(const RCSExportConfig& cfg);
	void mom_Dual_Pec_Cfie(const RCSExportConfig& cfg);
	void mom_Mono_Pec_Cfie(const RCSExportConfig& cfg);
	void mom_Dual_Die_Pmchwt(const RCSExportConfig& cfg);
	void mom_Mono_Die_Pmchwt(const RCSExportConfig& cfg);

	void exportZ(const RCSExportConfig& cfg);
};
#endif // MOM_H