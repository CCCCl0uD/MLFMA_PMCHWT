// ProcessFMM.h
#pragma once
#ifndef PROCESSFMM_H
#define PROCESSFMM_H

#include "GaussPoints.h"
#include "OverloadAlgo.h"
#include "OCTree.h"

#include <thread>
#include <atomic>
#include <iostream>
#include <omp.h>

namespace ProcessFMM {
	template<typename Kernel>
	inline std::vector<std::vector<Complex3D>> computeKernelMatrix(
		const std::vector<OCTree::Node*>& nodes, const std::vector<kp_Point>& kp_,
		const std::complex<double> k_wave, const int rows, const gaussPoints& gp,
		const double alpha, const std::complex<double> eta_,
		bool showProgressBar, Kernel&& kernelFunc)
	{
		std::vector<std::vector<Complex3D>> results(rows, std::vector<Complex3D>(kp_.size()));

		std::atomic<int> finishedRWG(0);
		std::thread progressThread;

		if (showProgressBar) {
			progressThread = std::thread([&]() {
				while (true) {
					int progressCount = finishedRWG.load(std::memory_order_relaxed);
					float progress = float(progressCount) / rows;
					std::cout << "\r[";
					int pos = int(50 * progress);
					for (int i = 0; i < 50; ++i) {
						if (i < pos) std::cout << "=";
						else if (i == pos) std::cout << ">";
						else std::cout << " ";
					}
					std::cout << "] " << std::setw(3) << int(progress * 100.0) << "%";
					std::cout.flush();
					if (progressCount >= rows) break;
					std::this_thread::sleep_for(std::chrono::milliseconds(200));
				}
				});
		}

#pragma omp parallel for schedule(dynamic)
		for (int nodeIdx = 0; nodeIdx < nodes.size(); ++nodeIdx) {
			OCTree::Node* node = nodes[nodeIdx];
			if (!node)continue;
			for (size_t rwgIdx = 0; rwgIdx < node->rwgIndices.size(); ++rwgIdx) {
				RWGBase* rwg = node->rwgIndices[rwgIdx];
				if (!rwg) continue;
				int id = rwg->rwgid;
				double l_ = rwg->edgeLength;
				double center_[3] = { node->center.x, node->center.y, node->center.z };
				const Point& freeP = *(rwg->freeVertexPositive);
				const Point& freeN = *(rwg->freeVertexNegative);
				const double normalP[3] = {
					rwg->positiveFace->externalNormalVector.x,
					rwg->positiveFace->externalNormalVector.y,
					rwg->positiveFace->externalNormalVector.z
				};
				const double normalN[3] = {
					rwg->negativeFace->externalNormalVector.x,
					rwg->negativeFace->externalNormalVector.y,
					rwg->negativeFace->externalNormalVector.z
				};
				const Point& vp0 = *(rwg->positiveFace->vertex1);
				const Point& vp1 = *(rwg->positiveFace->vertex2);
				const Point& vp2 = *(rwg->positiveFace->vertex3);
				const Point& vn0 = *(rwg->negativeFace->vertex1);
				const Point& vn1 = *(rwg->negativeFace->vertex2);
				const Point& vn2 = *(rwg->negativeFace->vertex3);

				for (size_t kpn = 0; kpn < kp_.size(); ++kpn) {
					double kp[3] = { kp_[kpn].x, kp_[kpn].y, kp_[kpn].z };
					Complex3D contrib{ 0.0,0.0,0.0 };
					for (int i = 0; i < gp.N_points_; ++i) {
						double w = gp.weight[i] * l_;

						double xpos = gp.l1[i] * vp0.x + gp.l2[i] * vp1.x + gp.l3[i] * vp2.x;
						double ypos = gp.l1[i] * vp0.y + gp.l2[i] * vp1.y + gp.l3[i] * vp2.y;
						double zpos = gp.l1[i] * vp0.z + gp.l2[i] * vp1.z + gp.l3[i] * vp2.z;
						Point rpos{ -1, xpos, ypos, zpos };

						double xneg = gp.l1[i] * vn0.x + gp.l2[i] * vn1.x + gp.l3[i] * vn2.x;
						double yneg = gp.l1[i] * vn0.y + gp.l2[i] * vn1.y + gp.l3[i] * vn2.y;
						double zneg = gp.l1[i] * vn0.z + gp.l2[i] * vn1.z + gp.l3[i] * vn2.z;
						Point rneg{ -1, xneg, yneg, zneg };

						Complex3D posContrib = kernelFunc(rpos, center_, freeP, kp, k_wave, w, normalP, alpha, eta_);
						Complex3D negContrib = kernelFunc(rneg, center_, freeN, kp, k_wave, w, normalN, alpha, eta_);

						contrib.x += posContrib.x - negContrib.x;
						contrib.y += posContrib.y - negContrib.y;
						contrib.z += posContrib.z - negContrib.z;
					}
					results[id][kpn] = contrib;
				}
				if (showProgressBar) {
					finishedRWG.fetch_add(1, std::memory_order_relaxed);
				}
			}
		}
		if (showProgressBar) {
			progressThread.join();
			std::cout << "\nInfo: Computation done, processed RWG bases: " << rows << " (Preprocessing.h)" << std::endl;
		}
		return results;
	}

	struct KernelVsmi_Efie {
		inline Complex3D operator()(const Point& r, const double center[3], const Point& freeVertex,
			double kp[3], const std::complex<double> k_wave, const double w, const double n_[3], const double alpha, const std::complex<double> eta_) const {
			// 向量 p = center - r
			double p[3] = { center[0] - r.x, center[1] - r.y, center[2] - r.z };
			double p_f[3] = { r.x - freeVertex.x, r.y - freeVertex.y, r.z - freeVertex.z };
			std::complex<double> phase = k_wave * (p[0] * kp[0] + p[1] * kp[1] + p[2] * kp[2]);
			std::complex<double> e = exp(-J * phase);
			Complex3D res;
			res.x = w * e * p_f[0];
			res.y = w * e * p_f[1];
			res.z = w * e * p_f[2];
			return res;
		}
	};

	struct KernelVfmj_Efie {
		inline Complex3D operator()(const Point& r, const double center[3], const Point& freeVertex,
			double kp[3], const std::complex<double> k_wave, const double w, const double n_[3], const double alpha, const std::complex<double> eta_) const {
			double p[3] = { r.x - center[0], r.y - center[1], r.z - center[2] };
			double p_f[3] = { r.x - freeVertex.x, r.y - freeVertex.y, r.z - freeVertex.z };
			std::complex<double> phase = k_wave * (p[0] * kp[0] + p[1] * kp[1] + p[2] * kp[2]);
			std::complex<double> e = std::exp(-J * phase);
			// 计算 - (k × (k × p_f))
			double temp[3], temp1[3];
			cross(temp, kp, p_f);
			cross(temp1, kp, temp);
			Complex3D res;
			res.x = w * e * (-temp1[0]);
			res.y = w * e * (-temp1[1]);
			res.z = w * e * (-temp1[2]);
			return res;
		}
	};

	struct KernelVfmj_Cfie {
		inline Complex3D operator()(const Point& r, const double center[3], const Point& freeVertex,
			double kp[3], const std::complex<double> k_wave, const double w, const double n_[3], const double alpha, const std::complex<double> eta_) const {
			double p[3] = { r.x - center[0], r.y - center[1], r.z - center[2] };
			double p_f[3] = { r.x - freeVertex.x, r.y - freeVertex.y, r.z - freeVertex.z };
			std::complex<double> phase = k_wave * (p[0] * kp[0] + p[1] * kp[1] + p[2] * kp[2]);
			std::complex<double> e = std::exp(-J * phase);
			// ---------- 第一项：-k × (k × p_f) ----------
			double k_cross_pf[3], k_cross_k_cross_pf[3];
			cross(k_cross_pf, kp, p_f);
			cross(k_cross_k_cross_pf, kp, k_cross_pf);
			double efie[3] = {
				-k_cross_k_cross_pf[0],
				-k_cross_k_cross_pf[1],
				-k_cross_k_cross_pf[2]
			};
			// ---------- 第二项：-k × (n × p_f) ----------
			double n_cross_pf[3], k_cross_n_cross_pf[3];
			cross(n_cross_pf, n_, p_f);
			cross(k_cross_n_cross_pf, kp, n_cross_pf);
			double mfie[3] = {
				k_cross_n_cross_pf[0],
				k_cross_n_cross_pf[1],
				k_cross_n_cross_pf[2]
			};

			Complex3D res;
			res.x = w * e * (alpha * efie[0] + (1 - alpha) * mfie[0]);
			res.y = w * e * (alpha * efie[1] + (1 - alpha) * mfie[1]);
			res.z = w * e * (alpha * efie[2] + (1 - alpha) * mfie[2]);
			return res;
		}
	};

	// 原始 FMM 中的 computeVsmi_Vfmi_Efie 可以改写为两个独立的调用
	inline void computeVsmi_Vfmj_Efie(
		std::vector<std::vector<Complex3D>>& Vsmi,
		std::vector<std::vector<Complex3D>>& Vfmj,
		const std::vector<OCTree::Node*>& nodes,
		const std::vector<kp_Point>& kp_, const EMSource& wave,
		const int rows, const gaussPoints& gausspoint, const double alpha,
		bool showProgress = true)
	{
		auto t_start = std::chrono::high_resolution_clock::now();
		Vsmi = computeKernelMatrix(nodes, kp_, wave.k1(), rows, gausspoint, alpha, wave.eta1(), showProgress, KernelVsmi_Efie());
		Vfmj = computeKernelMatrix(nodes, kp_, wave.k1(), rows, gausspoint, alpha, wave.eta1(), showProgress, KernelVfmj_Efie());
		auto t_end = std::chrono::high_resolution_clock::now();
		auto t_elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(t_end - t_start);
		std::cout << "Info: Vsmi/Vfmj (EFIE) done, elapsed: " << t_elapsed.count() / 1000.0
			<< " ms (ProcessFMM.h)" << std::endl;
	}

	inline void computeVsmi_Vfmj_Cfie(
		std::vector<std::vector<Complex3D>>& Vsmi,
		std::vector<std::vector<Complex3D>>& Vfmj,
		const std::vector<OCTree::Node*>& nodes,
		const std::vector<kp_Point>& kp_, const EMSource& wave,
		const int rows, const gaussPoints& gausspoint,
		const double alpha, bool showProgress = true)
	{
		auto t_start = std::chrono::high_resolution_clock::now();
		Vsmi = computeKernelMatrix(nodes, kp_, wave.k1(), rows, gausspoint, alpha, wave.eta1(), showProgress, KernelVsmi_Efie());
		Vfmj = computeKernelMatrix(nodes, kp_, wave.k1(), rows, gausspoint, alpha, wave.eta1(), showProgress, KernelVfmj_Cfie());
		auto t_end = std::chrono::high_resolution_clock::now();
		auto t_elapsed = std::chrono::duration_cast<std::chrono::microseconds>(t_end - t_start);
		std::cout << "Info: Vsmi/Vfmj (CFIE) done, elapsed: " << t_elapsed.count() / 1000.0
			<< " ms (ProcessFMM.h)" << std::endl;
	}

	inline void computeVsmi_Vfmj_Die(
		std::vector<std::vector<Complex3D>>& Vsmi,
		std::vector<std::vector<Complex3D>>& Vfmj,
		const std::vector<OCTree::Node*>& nodes,
		const std::vector<kp_Point>& kp_,
		const std::complex<double> k_wave, const std::complex<double> eta_,
		const int rows, const gaussPoints& gausspoint,
		const double alpha, bool showProgress = true)
	{
		auto t_start = std::chrono::high_resolution_clock::now();
		Vsmi = computeKernelMatrix(nodes, kp_, k_wave, rows, gausspoint, alpha, eta_, showProgress, KernelVsmi_Efie());
		Vfmj = computeKernelMatrix(nodes, kp_, k_wave, rows, gausspoint, alpha, eta_, showProgress, KernelVfmj_Efie());
		auto t_end = std::chrono::high_resolution_clock::now();
		auto t_elapsed = std::chrono::duration_cast<std::chrono::microseconds>(t_end - t_start);
		std::cout << "Info: Vsmi/Vfmj (Dielectric) done, elapsed: " << t_elapsed.count() / 1000.0
			<< " ms (ProcessFMM.h)" << std::endl;
	}
}

#endif // !PROCESSFMM_H
