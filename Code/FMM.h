// FMM.h
#pragma once
#ifndef FMM_H
#define FMM_H

#include "OCTree.h"
#include "RCSExportConfig.h"
#include "EMSource.h"
#include "GaussPoints.h"

struct VecHash {
	size_t operator()(const OCTree::nodePoint& v) const {
		return std::hash<double>()(v.x)
			^ std::hash<double>()(v.y)
			^ std::hash<double>()(v.z);
	}
};

class FMM {
public:
	FMM(const RCSExportConfig& cfg, const int selectIntegralEqu, const int selectMatrixSolver,
		std::string selectMono_Dual, const std::string pol_wave,
		const std::vector<std::vector<OCTree::Node*>>& octreeNodes,
		const std::vector<RWGBase>& rwgs, const int maxLevel_,
		const gaussPoints& gausspoint, const double E0, const EMSource& wave);

	const std::vector<std::vector<OCTree::Node*>>& octreeNodes_;    // Octree节点索引
	const std::vector<RWGBase>& rwgs;
	const int maxLevel_;                                            // 八叉树最大层级
	const double E0;
	const gaussPoints& gausspoint;
	const EMSource& wave;
	const int integralEquType_;
	const int matrixSolverType_;
	int L_k1;                                                // 多极子模式数
	int L_k2;
	int levelSpan;                                        // 级数跨度
	int row;                                              // RWG基函数个数
	/****************************************************************************************************************************/
	// === k1 ===
	std::vector<std::vector<kp_Point>> kp_lvl_k1;
	std::vector<std::vector<double>> theta_level_k1;         // 各层\theta取值 (maxLevel_, maxLevel_-1,..., 3, 2)
	std::vector<std::vector<double>> phi_level_k1;           // 各层\phi取值 (maxLevel_, maxLevel_-1,..., 3, 2)
	std::vector<double> WGL_phi_k1;                          // 存储各层的\phi方向权重
	std::vector<std::vector<double>> WGL_k1;                 // 按照高斯-勒让德积分点取点的各点权重

	// === k2 角度网格 ===
	std::vector<std::vector<kp_Point>> kp_lvl_k2;
	std::vector<std::vector<double>> theta_level_k2;         // 各层\theta取值 (maxLevel_, maxLevel_-1,..., 3, 2)
	std::vector<std::vector<double>> phi_level_k2;           // 各层\phi取值 (maxLevel_, maxLevel_-1,..., 3, 2)
	std::vector<double> WGL_phi_k2;                          // 存储各层的\phi方向权重
	std::vector<std::vector<double>> WGL_k2;                 // 按照高斯-勒让德积分点取点的各点权重
	/****************************************************************************************************************************/
	// === k1 ===
	std::vector<std::vector<Complex3D>> Vsmi, Vfmj;                             // 聚合因子矩阵的三维复数数组 (maxLevel_, maxLevel_-1,..., 3, 2)
	std::vector<std::vector<std::complex<double>>> TF_fmm;

	// === k2 ===
	std::vector<std::vector<Complex3D>> Vsmi2, Vfmj2;
	std::vector<std::vector<std::complex<double>>> TF_fmm2;

	std::vector<std::vector<OCTree::Node*>> node_far;
	std::vector<std::vector<int>> node_far_id, node_farVec_id;
	std::vector<OCTree::nodePoint> farVec;

	std::vector<std::vector<int>> Z_near_id;
	std::vector<std::vector<std::complex<double>>> Z_near;                      // 储存近场矩阵Z_near

	std::vector<std::vector<Complex3D>> Sm_fmm, Bm_fmm;

	std::vector<std::complex<double>> Vm;

	size_t computeMem();
	void matrix_solver(int n, std::complex<double> x[], std::complex<double> rhs[], int itr_max, int mr, double tol_abs, double tol_rel);

private:

	void fmm_Dual_Pec_Efie(const RCSExportConfig& cfg, const std::string pol_wave);
	void fmm_Mono_Pec_Efie(const RCSExportConfig& cfg, const std::string pol_wave);
	void fmm_Dual_Pec_Cfie(const RCSExportConfig& cfg, const std::string pol_wave);
	void fmm_Mono_Pec_Cfie(const RCSExportConfig& cfg, const std::string pol_wave);
	void fmm_Dual_Die_Pmchwt(const RCSExportConfig& cfg, const std::string pol_wave);
	void fmm_Mono_Die_Pmchwt(const RCSExportConfig& cfg, const std::string pol_wave);
};

#endif // FMM_H