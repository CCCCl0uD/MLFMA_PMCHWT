// ProcessMLFMM.h
#pragma once
#ifndef PROCESSMLFMM_H
#define PROCESSMLFMM_H

#include "OverloadAlgo.h"
#include "OCTree.h"

#include <omp.h>

namespace ProcessMLFMM {
	// ============================================================================
	// 1. clearSmBm — 清零多层 Sm/Bm 数组
	// ============================================================================
	inline void clearSmBm(std::vector<std::vector<std::vector<Complex3D>>>& Sm,
		std::vector<std::vector<std::vector<Complex3D>>>& Bm)
	{
#pragma omp parallel for
		for (int level = 0; level < static_cast<int>(Sm.size()); ++level) {
			int nodes = static_cast<int>(Sm[level].size());
			for (int node = 0; node < nodes; ++node) {
				int kp = static_cast<int>(Sm[level][node].size());
				for (int k = 0; k < kp; ++k) {
					Sm[level][node][k] = Complex3D{};
					Bm[level][node][k] = Complex3D{};
				}
			}
		}
	}

	// ============================================================================
	// 2. mexp — 最细层聚合: x × Vsmi → Sm[maxLevel_-2]
	// ============================================================================
	inline void mexp(std::vector<std::vector<std::vector<Complex3D>>>& Sm,
		const std::vector<OCTree::Node*>& finestNodes,
		const std::vector<std::vector<Complex3D>>& Vsmi,
		const std::vector<kp_Point>& kpFinest,
		const std::complex<double>* x, int xOffset,
		int maxLevel_)
	{
		const int nodes = static_cast<int>(finestNodes.size());
		const int maxlvl = maxLevel_ - 2;
		const int kps = static_cast<int>(kpFinest.size());
#pragma omp parallel for
		for (int nodeNum = 0; nodeNum < nodes; ++nodeNum) {
			OCTree::Node* node = finestNodes[nodeNum];
			int nRwg = static_cast<int>(node->rwgIndices.size());
			for (int ri = 0; ri < nRwg; ++ri) {
				int id = node->rwgIndices[ri]->rwgid;
				for (int ki = 0; ki < kps; ++ki) {
					Sm[maxlvl][nodeNum][ki].x += Vsmi[id][ki].x * x[xOffset + id];
					Sm[maxlvl][nodeNum][ki].y += Vsmi[id][ki].y * x[xOffset + id];
					Sm[maxlvl][nodeNum][ki].z += Vsmi[id][ki].z * x[xOffset + id];
				}
			}
		}
	}

	// ============================================================================
	// 3. m2m — 多极到多极 (子→父, 上行): 插值 + exp(-i k k̂·(r_father - r_child))
	// ============================================================================
	inline void m2m(std::vector<std::vector<std::vector<Complex3D>>>& Sm,
		const std::vector<std::vector<OCTree::Node*>>& octreeNodes,
		const std::vector<std::vector<kp_Point>>& kp_lvl,
		const std::vector<std::vector<InterpData>>& interpol,
		int maxLevel_, int levelSpan, std::complex<double> k_wave)
	{
		for (int level = maxLevel_ - 3; level > -1; --level) {
			const int rev = levelSpan - 1 - level;
			const int nodes = static_cast<int>(octreeNodes[level + 2].size());
			const int kps = static_cast<int>(kp_lvl[rev].size());
#pragma omp parallel for schedule(dynamic)
			for (int nodeNum = 0; nodeNum < nodes; ++nodeNum) {
				OCTree::Node* node = octreeNodes[level + 2][nodeNum];
				int sons = static_cast<int>(node->realSons.size());
				for (int sonN = 0; sonN < sons; ++sonN) {
					OCTree::Node* son = node->realSons[sonN];
					int id = son->nodeNumber;
					double vec2F[3] = {
						node->center.x - son->center.x,
						node->center.y - son->center.y,
						node->center.z - son->center.z
					};
					for (int kpn = 0; kpn < kps; ++kpn) {
						const InterpData& data = interpol[level][kpn];
						Complex3D f[3][3] = {
							{ Sm[level + 1][id][data.tp[0]], Sm[level + 1][id][data.tp[1]], Sm[level + 1][id][data.tp[2]] },
							{ Sm[level + 1][id][data.tp[3]], Sm[level + 1][id][data.tp[4]], Sm[level + 1][id][data.tp[5]] },
							{ Sm[level + 1][id][data.tp[6]], Sm[level + 1][id][data.tp[7]], Sm[level + 1][id][data.tp[8]] }
						};
						const auto& kp_ = kp_lvl[rev][kpn];
						std::complex<double> kr = k_wave * (vec2F[0] * kp_.x + vec2F[1] * kp_.y + vec2F[2] * kp_.z);
						std::complex<double> e_ = exp(-J * kr);
						Complex3D result{};
						for (int tt = 0; tt < 3; ++tt)
							for (int pp = 0; pp < 3; ++pp) {
								result.x += f[tt][pp].x * data.w[tt * 3 + pp];
								result.y += f[tt][pp].y * data.w[tt * 3 + pp];
								result.z += f[tt][pp].z * data.w[tt * 3 + pp];
							}
						Sm[level][nodeNum][kpn].x += e_ * result.x;
						Sm[level][nodeNum][kpn].y += e_ * result.y;
						Sm[level][nodeNum][kpn].z += e_ * result.z;
					}
				}
			}
		}
	}

	// ============================================================================
	// 4. m2l — 多极到局部 (转移): Bm = w · T · Sm  (Sm 只读, 权重在此施加)
	// ============================================================================
	inline void m2l(const std::vector<std::vector<std::vector<Complex3D>>>& Sm,
		std::vector<std::vector<std::vector<Complex3D>>>& Bm,
		const std::vector<std::vector<OCTree::Node*>>& octreeNodes,
		const std::vector<std::vector<std::vector<std::complex<double>>>>& TF,
		const std::vector<std::vector<kp_Point>>& kp_lvl,
		const std::vector<std::vector<double>>& WGL,
		const std::vector<double>& WGL_phi,
		const std::vector<std::vector<double>>& phi_level,
		int maxLevel_, int levelSpan)
	{
		for (int level = 0; level < maxLevel_ - 1; ++level) {
			const int rev = levelSpan - 1 - level;
			const int kps = static_cast<int>(kp_lvl[rev].size());
			const int nodes = static_cast<int>(octreeNodes[level + 2].size());
			const int pN = static_cast<int>(phi_level[rev].size());
#pragma omp parallel for
			for (int nodeNum = 0; nodeNum < nodes; ++nodeNum) {
				OCTree::Node* node = octreeNodes[level + 2][nodeNum];
				const int drN = static_cast<int>(node->distantRelatives.size());
				auto& Bmvec = Bm[level][nodeNum];
				for (int drIdx = 0; drIdx < drN; ++drIdx) {
					int vecid = node->dRVecId[drIdx];
					int drid = node->distantRelatives[drIdx]->nodeNumber;
					auto& TFvec = TF[level][vecid];
					auto& Smvec = Sm[level][drid];
					for (int kp_i = 0; kp_i < kps; ++kp_i) {
						int t = kp_i / pN;
						double weight = WGL[rev][t] * WGL_phi[rev];
						Bmvec[kp_i].x += TFvec[kp_i] * Smvec[kp_i].x * weight;
						Bmvec[kp_i].y += TFvec[kp_i] * Smvec[kp_i].y * weight;
						Bmvec[kp_i].z += TFvec[kp_i] * Smvec[kp_i].z * weight;
					}
				}
			}
		}
	}

	// ============================================================================
	// 5. l2l — 局部到局部 (父→子, 下行): exp(-i k k̂·(r_child - r_father)) + 反插值
	// ============================================================================
	inline void l2l(std::vector<std::vector<std::vector<Complex3D>>>& Bm,
		const std::vector<std::vector<OCTree::Node*>>& octreeNodes,
		const std::vector<std::vector<kp_Point>>& kp_lvl,
		const std::vector<std::vector<double>>& phi_level,
		const std::vector<std::vector<InterpData>>& interpol,
		int maxLevel_, int levelSpan, std::complex<double> k_wave)
	{
		for (int level = 0; level < maxLevel_ - 2; ++level) {
			const int revF = levelSpan - 1 - level;
			const int revS = levelSpan - 2 - level;
			auto& nodesF = octreeNodes[level + 2];
			auto& kpLevel = kp_lvl[revF];
			auto& interpLevel = interpol[level];
			const int nodefN = static_cast<int>(nodesF.size());
			const int kpf = static_cast<int>(kpLevel.size());
			const int phiC = static_cast<int>(phi_level[revS].size());
#pragma omp parallel for
			for (int fIdx = 0; fIdx < nodefN; ++fIdx) {
				OCTree::Node* nodef = nodesF[fIdx];
				int sonN = static_cast<int>(nodef->realSons.size());
				for (int son = 0; son < sonN; ++son) {
					double r[3]{
						nodef->realSons[son]->center.x - nodef->center.x,
						nodef->realSons[son]->center.y - nodef->center.y,
						nodef->realSons[son]->center.z - nodef->center.z
					};
					int sonid = nodef->realSons[son]->nodeNumber;
					auto& Bm1Vec = Bm[level + 1][sonid];
					for (int kp_i = 0; kp_i < kpf; ++kp_i) {
						auto& kp_ = kpLevel[kp_i];
						std::complex<double> kr = k_wave * (kp_.x * r[0] + kp_.y * r[1] + kp_.z * r[2]);
						std::complex<double> e_ = exp(-J * kr);
						const InterpData& data = interpLevel[kp_i];
						auto& Bm2Vec = Bm[level][fIdx][kp_i];
						for (int tt = 0; tt < 3; ++tt) {
							int tIndex = data.tIdx[tt];
							for (int pp = 0; pp < 3; ++pp) {
								int pIndex = data.pIdx[pp];
								int index = tIndex * phiC + pIndex;
								auto& dst = Bm1Vec[index];
								double wgt = data.w[tt * 3 + pp];
								dst.x += e_ * Bm2Vec.x * wgt;
								dst.y += e_ * Bm2Vec.y * wgt;
								dst.z += e_ * Bm2Vec.z * wgt;
							}
						}
					}
				}
			}
		}
	}

	// ============================================================================
	// 6. lexpKpDot — 单个 kp 方向 Bm · Vfmj 内积
	// ============================================================================
	inline std::complex<double> lexpKpDot(const Complex3D& B, const Complex3D& V) {
		return V.x * B.x + V.y * B.y + V.z * B.z;
	}

	// ============================================================================
	// 7. lexpPEC — PEC 解聚合 (最细层) + 近场
	// ============================================================================
	inline void lexpPEC(std::complex<double>* w,
		const std::vector<OCTree::Node*>& finestNodes,
		const std::vector<std::vector<Complex3D>>& Vfmj,
		const std::vector<std::vector<std::vector<Complex3D>>>& Bm,
		int tN, int pN, int maxLevel_,
		std::complex<double> k_, std::complex<double> eta_,
		const std::vector<std::vector<std::complex<double>>>& Z_near,
		const std::vector<std::vector<int>>& Z_near_id,
		const std::complex<double>* x)
	{
		const int nodeNum = static_cast<int>(finestNodes.size());
#pragma omp parallel for schedule(dynamic)
		for (int nodeIdx = 0; nodeIdx < nodeNum; ++nodeIdx) {
			OCTree::Node* node = finestNodes[nodeIdx];
			const int rwgNum = static_cast<int>(node->rwgIndices.size());
			const auto& B = Bm[maxLevel_ - 2][nodeIdx];
			for (int r = 0; r < rwgNum; ++r) {
				int rwgid = node->rwgIndices[r]->rwgid;
				const auto& V = Vfmj[rwgid];
				std::complex<double> val(0.0, 0.0);
				for (int tt = 0; tt < tN; ++tt)
					for (int pp = 0; pp < pN; ++pp)
						val += lexpKpDot(B[tt * pN + pp], V[tt * pN + pp]);
				val *= (k_ * eta_ / (4.0 * Pi)) * (k_ / (4.0 * Pi));
				const int near = static_cast<int>(Z_near[rwgid].size());
				for (int j = 0; j < near; ++j)
					val += x[Z_near_id[rwgid][j]] * Z_near[rwgid][j];
				w[rwgid] = val;
			}
		}
	}

	// ============================================================================
	// 8. lexpPMCHWT — PMCHWT 单次 Pass 解聚合 (L + K 同时计算)
	//    constL/constK : L/K 算子常数 (含符号)
	//    wOffset_L/wOffset_K : 写入 w 的偏移 (0 或 N)
	// ============================================================================
	inline void lexpPMCHWT(std::complex<double>* w,
		std::complex<double> constL, std::complex<double> constK,
		int wOffset_L, int wOffset_K, int N,
		const std::vector<OCTree::Node*>& finestNodes,
		const std::vector<std::vector<Complex3D>>& Vfmj,
		const std::vector<std::vector<std::vector<Complex3D>>>& Bm,
		const std::vector<kp_Point>& kpFinest,
		int tN, int pN, int maxLevel_)
	{
		const int nodeNum = static_cast<int>(finestNodes.size());
		const int finestIdx = maxLevel_ - 2;
#pragma omp parallel for schedule(dynamic)
		for (int nodeIdx = 0; nodeIdx < nodeNum; ++nodeIdx) {
			OCTree::Node* nd = finestNodes[nodeIdx];
			const auto& BB = Bm[finestIdx][nodeIdx];
			for (int ri = 0; ri < static_cast<int>(nd->rwgIndices.size()); ++ri) {
				int rwgid = nd->rwgIndices[ri]->rwgid;
				const auto& V = Vfmj[rwgid];
				std::complex<double> valL(0.0, 0.0);
				std::complex<double> valK(0.0, 0.0);
				for (int tt = 0; tt < tN; ++tt) {
					for (int pp = 0; pp < pN; ++pp) {
						int idx = tt * pN + pp;
						valL += lexpKpDot(BB[idx], V[idx]);
						double kpp[3] = { kpFinest[idx].x, kpFinest[idx].y, kpFinest[idx].z };
						Complex3D kcv;
						cross(kcv, kpp, V[idx]);
						valK += lexpKpDot(BB[idx], kcv);
					}
				}
				w[wOffset_L + rwgid] += constL * valL;
				w[wOffset_K + rwgid] += constK * valK;
			}
		}
	}
} // namespace ProcessMLFMM
#endif // !PROCESSMLFMM_H
