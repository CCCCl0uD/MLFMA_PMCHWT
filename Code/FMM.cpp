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
		AlgoFMM::computeFMM_far(octreeNodes_[maxLevel_],
			farVec, node_far, node_far_id, node_farVec_id,
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
		// k1: TF_fmm1
		AlgoFMM::computeFMM_far(octreeNodes_[maxLevel_],
			farVec, node_far, node_far_id, node_farVec_id,
			TF_fmm, kp_lvl_k1, wave.k1(), L_k1);

		// k2: TF_fmm2
		AlgoFMM::computeFMM_far(octreeNodes_[maxLevel_],
			farVec, node_far, node_far_id, node_farVec_id,
			TF_fmm2, kp_lvl_k2, wave.k2(), L_k2);

		// k1: Vsmi1/Vfmj1
		ProcessFMM::computeVsmi_Vfmj_Die(Vsmi, Vfmj,
			octreeNodes_[maxLevel_], kp_lvl_k1[0], wave.k1(), wave.eta1(),
			row, gausspoint, alpha, 1);

		// k2: Vsmi2/Vfmj2
		ProcessFMM::computeVsmi_Vfmj_Die(Vsmi2, Vfmj2,
			octreeNodes_[maxLevel_], kp_lvl_k2[0], wave.k2(), wave.eta2(),
			row, gausspoint, alpha, 1);

		// PMCHWT: Z_near 2N×2N
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
		int K1 = static_cast<int>(kp_lvl_k1[0].size());

		if (EquType_ != 2) {
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
			// ===== PMCHW: 8 次子运算 =====
			const int N = row;
			const int K2 = static_cast<int>(kp_lvl_k2[0].size());
			const int tN1 = static_cast<int>(theta_level_k1[0].size());
			const int pN1 = static_cast<int>(phi_level_k1[0].size());
			const int tN2 = static_cast<int>(theta_level_k2[0].size());
			const int pN2 = static_cast<int>(phi_level_k2[0].size());
			const std::complex<double> const1 = wave.k1() * wave.k1() * wave.eta1() / (16.0 * Pi * Pi);
			const std::complex<double> const2 = wave.k2() * wave.k2() / (16.0 * Pi * Pi);

			std::vector<std::vector<Complex3D>> Sm_J(nodenum, std::vector<Complex3D>(K1, Complex3D()));
			std::vector<std::vector<Complex3D>> Bm_J(nodenum, std::vector<Complex3D>(K1, Complex3D()));
			std::vector<std::vector<Complex3D>> Sm_M(nodenum, std::vector<Complex3D>(K2, Complex3D()));
			std::vector<std::vector<Complex3D>> Bm_M(nodenum, std::vector<Complex3D>(K2, Complex3D()));

			// ---- 区域1: Mexp(J) + ML(J) + Mexp(M) + ML(M) ----
#pragma omp parallel for schedule(dynamic)
			for (int nodeN = 0; nodeN < nodenum; ++nodeN) {
				OCTree::Node* nd = octreeNodes_[maxLevel_][nodeN];

				for (int rwg_i = 0; rwg_i < nd->rwgIndices.size(); ++rwg_i) {
					int id = nd->rwgIndices[rwg_i]->rwgid;

					for (int kpi = 0; kpi < K1; ++kpi) {
						Sm_J[nodeN][kpi].x += Vsmi[id][kpi].x * x[id];       // x[0..N-1]=I_J
						Sm_J[nodeN][kpi].y += Vsmi[id][kpi].y * x[id];
						Sm_J[nodeN][kpi].z += Vsmi[id][kpi].z * x[id];
						Sm_M[nodeN][kpi].x += Vsmi[id][kpi].x * x[N + id];  // x[N..2N-1]=I_M
						Sm_M[nodeN][kpi].y += Vsmi[id][kpi].y * x[N + id];
						Sm_M[nodeN][kpi].z += Vsmi[id][kpi].z * x[N + id];
					}
				}
			}

#pragma omp parallel for schedule(dynamic)
			for (int nodeN = 0; nodeN < nodenum; ++nodeN) {
				for (int farn = 0; farn < static_cast<int>(node_far[nodeN].size()); ++farn) {
					int farId = node_far_id[nodeN][farn];
					const auto& TF = TF_fmm[node_farVec_id[nodeN][farn]];

					for (int kpi = 0; kpi < K1; ++kpi) {
						Bm_J[nodeN][kpi].x += TF[kpi] * Sm_J[farId][kpi].x;
						Bm_J[nodeN][kpi].y += TF[kpi] * Sm_J[farId][kpi].y;
						Bm_J[nodeN][kpi].z += TF[kpi] * Sm_J[farId][kpi].z;
						Bm_M[nodeN][kpi].x += TF[kpi] * Sm_M[farId][kpi].x;
						Bm_M[nodeN][kpi].y += TF[kpi] * Sm_M[farId][kpi].y;
						Bm_M[nodeN][kpi].z += TF[kpi] * Sm_M[farId][kpi].z;
					}
				}
			}

			// ---- 区域1: Lexp (4 个子运算) + 区域2: Mexp/ML/Lexp ----
			// 为区域2初始化临时数组
			std::vector<std::vector<Complex3D>> Sm_J2(nodenum, std::vector<Complex3D>(K1, Complex3D()));
			std::vector<std::vector<Complex3D>> Bm_J2(nodenum, std::vector<Complex3D>(K1, Complex3D()));
			std::vector<std::vector<Complex3D>> Sm_M2(nodenum, std::vector<Complex3D>(K1, Complex3D()));
			std::vector<std::vector<Complex3D>> Bm_M2(nodenum, std::vector<Complex3D>(K1, Complex3D()));

			// Lexp1 (区域1 Lexp 需用 Bm_J 和 Bm_M):
			// 先完成解聚① L1(J), ② K1(J), ③ K1(M), ④ L1(M)
			// 同时做区域2的 Mexp(J), Mexp(M)
#pragma omp parallel for schedule(dynamic)
			for (int nodeN = 0; nodeN < nodenum; ++nodeN) {
				OCTree::Node* nd = octreeNodes_[maxLevel_][nodeN];
				for (int r = 0; r < static_cast<int>(nd->rwgIndices.size()); ++r) {
					int rwgid = nd->rwgIndices[r]->rwgid;

					// 区域2 Mexp(J)
					for (int kpi = 0; kpi < K2; ++kpi) {
						Sm_J2[nodeN][kpi].x += Vsmi2[rwgid][kpi].x * x[rwgid];
						Sm_J2[nodeN][kpi].y += Vsmi2[rwgid][kpi].y * x[rwgid];
						Sm_J2[nodeN][kpi].z += Vsmi2[rwgid][kpi].z * x[rwgid];
						Sm_M2[nodeN][kpi].x += Vsmi2[rwgid][kpi].x * x[N + rwgid];
						Sm_M2[nodeN][kpi].y += Vsmi2[rwgid][kpi].y * x[N + rwgid];
						Sm_M2[nodeN][kpi].z += Vsmi2[rwgid][kpi].z * x[N + rwgid];
					}

					// ---- ① L1(J): EFIE on J → w[0:N-1] ----
					std::complex<double> valL1J(0.0, 0.0);
					const auto& V1 = Vfmj[rwgid];
					for (int ti = 0; ti < theta_level_k1[0].size(); ++ti) {
						for (int pj = 0; pj < pN1; ++pj) {
							int idx = ti * pN1 + pj;
							double d2k = WGL_k1[0][ti] * WGL_phi_k1[0];
							valL1J += d2k * dot(V1[idx], Bm_J[nodeN][idx]);
						}
					}
					w[rwgid] += const1 * valL1J;

					// ---- ② K1(J): K on J → w[N:2N-1]（负号）----
					// Vfmj_K = k̂ × Vfmj（在 xyz 中直接用叉乘）
					std::complex<double> valK1J(0.0, 0.0);
					for (int ti = 0; ti < tN1; ++ti) {
						for (int pj = 0; pj < pN1; ++pj) {
							int idx = ti * pN1 + pj;
							double d2k = WGL_k1[0][ti] * WGL_phi_k1[0];
							double kp[3] = { kp_lvl_k1[0][idx].x, kp_lvl_k1[0][idx].y, kp_lvl_k1[0][idx].z };
							Complex3D kCrossVfmj;
							cross(kCrossVfmj, kp, V1[idx]);
							valK1J += d2k * dot(kCrossVfmj, Bm_J[nodeN][idx]);
						}
					}
					w[N + rwgid] -= const1 * valK1J;

					// ---- ③ K1(M): K on M → w[0:N-1] ----
					std::complex<double> valK1M(0.0, 0.0);
					for (int ti = 0; ti < tN1; ++ti) {
						for (int pj = 0; pj < pN1; ++pj) {
							int idx = ti * pN1 + pj;
							double d2k = WGL_k1[0][ti] * WGL_phi_k1[0];
							double kp[3] = { kp_lvl_k1[0][idx].x, kp_lvl_k1[0][idx].y, kp_lvl_k1[0][idx].z };
							Complex3D kCrossVfmj;
							cross(kCrossVfmj, kp, Vfmj[rwgid][idx]);
							valK1M += d2k * dot(kCrossVfmj, Bm_M[nodeN][idx]);
						}
					}
					w[rwgid] += const1 * valK1M;

					// ---- ④ L1(M): EFIE on M → w[N:2N-1] ----
					std::complex<double> valL1M(0.0, 0.0);
					for (int ti = 0; ti < tN1; ++ti) {
						for (int pj = 0; pj < pN1; ++pj) {
							int idx = ti * pN1 + pj;
							double d2k = WGL_k1[0][ti] * WGL_phi_k1[0];
							valL1M += d2k * dot(Vfmj[rwgid][idx], Bm_M[nodeN][idx]);
						}
					}
					w[N + rwgid] += const1 * valL1M;
				}
			}

			// ---- 区域2: ML(J) 和 ML(M) ----
#pragma omp parallel for schedule(dynamic)
			for (int nodeN = 0; nodeN < nodenum; ++nodeN) {
				for (int farn = 0; farn < node_far[nodeN].size(); ++farn) {
					int farId = node_far_id[nodeN][farn];
					const auto& TF2 = TF_fmm2[node_farVec_id[nodeN][farn]];
					for (int kpi = 0; kpi < K2; ++kpi) {
						Bm_J2[nodeN][kpi].x += TF2[kpi] * Sm_J2[farId][kpi].x;
						Bm_J2[nodeN][kpi].y += TF2[kpi] * Sm_J2[farId][kpi].y;
						Bm_J2[nodeN][kpi].z += TF2[kpi] * Sm_J2[farId][kpi].z;
						Bm_M2[nodeN][kpi].x += TF2[kpi] * Sm_M2[farId][kpi].x;
						Bm_M2[nodeN][kpi].y += TF2[kpi] * Sm_M2[farId][kpi].y;
						Bm_M2[nodeN][kpi].z += TF2[kpi] * Sm_M2[farId][kpi].z;
					}
				}
			}

			// ---- 区域2: Lexp (⑤~⑧) + 近场 ----
#pragma omp parallel for schedule(dynamic)
			for (int nodeIdx = 0; nodeIdx < nodenum; ++nodeIdx) {
				OCTree::Node* nd = octreeNodes_[maxLevel_][nodeIdx];
				for (int r = 0; r < nd->rwgIndices.size(); ++r) {
					int rwgid = nd->rwgIndices[r]->rwgid;
					const auto& V2 = Vfmj2[rwgid];
					// ⑤ L2(J): const2 * eta2 * L2(J) → w[0:N-1]
					std::complex<double> vL2J(0.0, 0.0);
					for (int ti = 0; ti < tN2; ++ti) {
						for (int pj = 0; pj < pN2; ++pj) {
							int idx = ti * pN2 + pj;
							double d2k = WGL_k2[0][ti] * WGL_phi_k2[0];
							vL2J += d2k * dot(V2[idx], Bm_J2[nodeIdx][idx]);
						}
					}
					w[rwgid] += const2 * wave.eta2() * vL2J;

					// ⑥ K2(J): -const2 * eta1 * K2(J) → w[N:2N-1]
					std::complex<double> vK2J(0.0, 0.0);
					for (int ti = 0; ti < tN2; ++ti) {
						for (int pj = 0; pj < pN2; ++pj) {
							int idx = ti * pN2 + pj;
							double d2k = WGL_k2[0][ti] * WGL_phi_k2[0];
							double kp[3] = { kp_lvl_k2[0][idx].x, kp_lvl_k2[0][idx].y, kp_lvl_k2[0][idx].z };
							Complex3D kcv;
							cross(kcv, kp, Vfmj2[rwgid][idx]);
							vK2J += d2k * dot(kcv, Bm_J2[nodeIdx][idx]);
						}
					}
					w[N + rwgid] -= const2 * wave.eta1() * vK2J;

					// ⑦ K2(M): const2 * eta1 * K2(M) → w[0:N-1]
					std::complex<double> vK2M(0.0, 0.0);
					for (int ti = 0; ti < tN2; ++ti) {
						for (int pj = 0; pj < pN2; ++pj) {
							int idx = ti * pN2 + pj;
							double d2k = WGL_k2[0][ti] * WGL_phi_k2[0];
							double kp[3] = { kp_lvl_k2[0][idx].x, kp_lvl_k2[0][idx].y, kp_lvl_k2[0][idx].z };
							Complex3D kcv;
							cross(kcv, kp, Vfmj2[rwgid][idx]);
							vK2M += d2k * dot(kcv, Bm_M2[nodeIdx][idx]);
						}
					}
					w[rwgid] += const2 * wave.eta1() * vK2M;

					// ⑧ L2(M): const2 * eta2 * εr * L2(M) → w[N:2N-1]
					std::complex<double> vL2M(0.0, 0.0);
					for (int ti = 0; ti < tN2; ++ti) {
						for (int pj = 0; pj < pN2; ++pj) {
							int idx = ti * pN2 + pj;
							double d2k = WGL_k2[0][ti] * WGL_phi_k2[0];
							vL2M += d2k * dot(Vfmj2[rwgid][idx], Bm_M2[nodeIdx][idx]);
						}
					}
					w[N + rwgid] += const2 * wave.eta2() * wave.epsilonR() * vL2M;

					// ---- 近场: 2N×2N 稀疏矩阵 × x ----
					// J 测试行
					for (int j = 0; j < static_cast<int>(Z_near[rwgid].size()); ++j) {
						w[rwgid] += x[Z_near_id[rwgid][j]] * Z_near[rwgid][j];
					}
					// M 测试行
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
			solver.wave.k1(), solver.wave.k2(), kInc, eInc, hInc, solver.Vm);
		};
	RCSUtils::computeDualStatic_PMCHWT(*this, cfg, computeV);
}

void FMM::fmm_Mono_Die_Pmchwt(const RCSExportConfig& cfg, const std::string pol_wave)
{
	auto computeV = [](FMM& solver, double kInc[3], double eInc[3], double hInc[3]) {
		RHS::computeV_PMCHWT(solver.rwgs, solver.gausspoint,
			solver.wave.k1(), solver.wave.k2(), kInc, eInc, hInc, solver.Vm);
		};
	RCSUtils::computeMonoStatic_PMCHWT(*this, cfg, computeV);
}