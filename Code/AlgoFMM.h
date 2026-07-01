// AlgoFMM.h
#pragma once
#ifndef ALGOFMM_H
#define ALGOFMM_H

#include <atomic>
#include <chrono>
#include <complex>
#include <iomanip>
#include <iostream>
#include <thread>
#include <unordered_set>
#include <vector>
#include <boost/math/special_functions/legendre.hpp>

#include "Bessel.h"
#include "OverloadAlgo.h"
#include "EMSource.h"
#include "OCTree.h"

namespace AlgoFMM {
	// 牛顿-拉夫逊方法来寻找勒让德多项式的根
	inline double newton_raphson(int n, double x0, double tolerance, int max_iter) {
		double x = x0;
		for (int i = 0; i < max_iter; ++i) {
			double p = boost::math::legendre_p(n, x);
			double dp = boost::math::legendre_p_prime(n, x);
			double dx = p / dp;
			x -= dx;
			if (std::abs(dx) < tolerance) {
				break;
			}
		}
		return x;
	}

	// 计算n点高斯-勒让德积分的节点和权重
	inline void gauss_legendre(int n, std::vector<double>& XGL_temp, std::vector<double>& WGL_temp) {
		XGL_temp.resize(n);
		WGL_temp.resize(n);
		const double tolerance = 1e-10;
		const int max_iter = 100;
		for (int i = 0; i < (n + 1) / 2; ++i) {
			// 初始猜测值
			double x0 = std::cos(Pi * (i + 0.75) / (n + 0.5));
			double root = newton_raphson(n, x0, tolerance, max_iter);
			XGL_temp[i] = -root;
			XGL_temp[n - i - 1] = root;

			double dp = boost::math::legendre_p_prime(n, root);
			double weight = 2.0 / ((1 - root * root) * dp * dp);
			WGL_temp[i] = weight;
			WGL_temp[n - i - 1] = weight;
		}
	}

	// FMM远场数据处理
	inline void computeFMMFarGroups(
		const std::vector<OCTree::Node*>& nodes,
		std::vector<OCTree::nodePoint>& farVec,
		std::vector<std::vector<OCTree::Node*>>& node_far,
		std::vector<std::vector<int>>& node_far_id,
		std::vector<std::vector<int>>& node_farVec_id)
	{
		farVec.clear();

		const int num = static_cast<int>(nodes.size());
		node_far.resize(num);
		node_far_id.resize(num);
		node_farVec_id.resize(num);

		for (int i = 0; i < num; ++i) {
			node_far[i].clear();
			node_far_id[i].clear();
			node_farVec_id[i].clear();
		}

		for (int i = 0; i < num; ++i) {
			OCTree::Node* node1 = nodes[i];

			std::unordered_set<OCTree::Node*> nearSet(
				node1->nearNeighbors.begin(),
				node1->nearNeighbors.end());
			nearSet.insert(node1);

			for (int j = 0; j < num; ++j) {
				OCTree::Node* node2 = nodes[j];
				if (nearSet.find(node2) == nearSet.end()) {
					node_far[i].push_back(node2);
					node_far_id[i].push_back(j);
				}
			}
		}

		for (int idx1 = 0; idx1 < num; ++idx1) {
			OCTree::Node* node = nodes[idx1];
			if (node_far[idx1].empty()) {
				continue;
			}

			for (int idx2 = 0; idx2 < static_cast<int>(node_far[idx1].size()); ++idx2) {
				OCTree::Node* far = node_far[idx1][idx2];

				OCTree::nodePoint vecToDR = {
					OCTree::snapZero(node->center.x - far->center.x),
					OCTree::snapZero(node->center.y - far->center.y),
					OCTree::snapZero(node->center.z - far->center.z)
				};

				int vecIdx = -1;
				for (int l = 0; l < static_cast<int>(farVec.size()); ++l) {
					if (OCTree::isVecEqual(vecToDR, farVec[l])) {
						vecIdx = l;
						break;
					}
				}

				if (vecIdx == -1) {
					vecIdx = static_cast<int>(farVec.size());
					farVec.push_back(vecToDR);
				}

				node_farVec_id[idx1].push_back(vecIdx);
			}
		}
	}

	inline void computeFMM_far(std::vector<OCTree::nodePoint>& farVec,
		std::vector<std::vector<OCTree::Node*>>& node_far,
		std::vector<std::vector<int>>& node_far_id,
		std::vector<std::vector<int>>& node_farVec_id,
		std::vector<std::vector<std::complex<double>>>& TF_fmm,
		std::vector<std::vector<kp_Point>>& kp_lvl,
		const std::complex<double> k_, const int L)
	{
		auto t_start = std::chrono::high_resolution_clock::now();

		// 初始化转移因子矩阵
		TF_fmm.resize(farVec.size());
		for (int i = 0; i < farVec.size(); ++i) {
			TF_fmm[i].assign(kp_lvl[0].size(), std::complex<double>(0.0, 0.0));
		}

		// 计算转移因子
		std::atomic<int> finishedNodes(0);
		const int totalNodes = farVec.size();
		const int barWidth = 50;

		std::thread progressThread([&]() {
			while (true) {
				int progressCount = finishedNodes.load(std::memory_order_relaxed);
				float progress = float(progressCount) / totalNodes;

				std::cout << "\r[";
				int pos = int(barWidth * progress);

				for (int i = 0; i < barWidth; ++i) {
					if (i < pos) std::cout << "=";
					else if (i == pos) std::cout << ">";
					else std::cout << " ";
				}

				std::cout << "] "
					<< std::setw(3)
					<< int(progress * 100.0) << "%";

				std::cout.flush();

				if (progressCount >= totalNodes)
					break;

				std::this_thread::sleep_for(std::chrono::milliseconds(300));
			}
			});

#pragma omp parallel
		{
#pragma omp for schedule(dynamic)
			for (int nodeFarIdx = 0; nodeFarIdx < farVec.size(); ++nodeFarIdx) {
				double r[3] = {
					farVec[nodeFarIdx].x,
					farVec[nodeFarIdx].y,
					farVec[nodeFarIdx].z
				};
				double Rpq = std::sqrt(r[0] * r[0] + r[1] * r[1] + r[2] * r[2]);
				normalize(r, Rpq);

				std::complex<double> kR = k_ * Rpq;

				const int kps = kp_lvl[0].size();

				for (int kpNum = 0; kpNum < kps; ++kpNum) {
					double kp[3] = {
						kp_lvl[0][kpNum].x,
						kp_lvl[0][kpNum].y,
						kp_lvl[0][kpNum].z
					};
					// 勒让德多项式的参数
					double r_dot_kp = dot(r, kp);
					std::vector<std::complex<double>> h2_array;
					ComplexBessel::spherical_hankel_2_array(L, kR, h2_array);

					for (int l = 0; l <= L; ++l) {
						double p = boost::math::legendre_p(l, r_dot_kp);
						std::complex<double> jl = jl_table[l & 3];
						TF_fmm[nodeFarIdx][kpNum] += jl * h2_array[l] * (2.0 * l + 1.0) * p;
					}
				}
				finishedNodes.fetch_add(1, std::memory_order_relaxed);
			}
		}

		progressThread.join();
		auto t_end = std::chrono::high_resolution_clock::now();
		auto t_elapsed = std::chrono::duration_cast<std::chrono::microseconds>(t_end - t_start);
		std::cout << "\nInfo: Transfer factor calculation done, elapsed: "
			<< t_elapsed.count() / 1000.0 << " ms (AlgoFMM.h)\n";
	}

	inline void computeBase(
		const std::vector<std::vector<OCTree::Node*>>& octreeNodes,
		const EMSource& wave, const int integralEquType_,
		int maxLevel_, int rwgSize, int& L_k1, int& L_k2, int& row, int& levelSpan,
		std::vector<std::vector<double>>& WGL_k1, std::vector<std::vector<double>>& WGL_k2,
		std::vector<double>& WGL_phi_k1, std::vector<double>& WGL_phi_k2,
		std::vector<std::vector<double>>& theta_level_k1, std::vector<std::vector<double>>& theta_level_k2,
		std::vector<std::vector<double>>& phi_level_k1, std::vector<std::vector<double>>& phi_level_k2,
		std::vector<std::vector<kp_Point>>& kp_lvl_k1, std::vector<std::vector<kp_Point>>& kp_lvl_k2)
	{
		const double highLevelLen = octreeNodes[maxLevel_][0]->cubeLength;

		double D_ = std::sqrt(3.0) * highLevelLen;
		L_k1 = static_cast<int>(std::ceil(wave.k1_abs() * D_ + 6.0 * std::cbrt(wave.k1_abs() * D_)));

		row = rwgSize;
		levelSpan = maxLevel_ - 1;

		// 局部临时变量
		std::vector<double> XGL_temp, WGL_temp;
		std::vector<double> theta, phi;
		std::vector<std::vector<double>> XGL_local;

		for (int i = 0; i < levelSpan; i++) {
			int n = (1 << i) * L_k1;
			AlgoFMM::gauss_legendre(n, XGL_temp, WGL_temp);
			XGL_local.push_back(XGL_temp);
			WGL_k1.push_back(WGL_temp);
			XGL_temp.clear();
			WGL_temp.clear();
		}

		for (int i = 0; i < levelSpan; i++) {
			theta.clear();
			phi.clear();

			for (int j = static_cast<int>(WGL_k1[i].size()) - 1; j >= 0; j--) {
				double xGL = XGL_local[i][j];  // 原来从 XGL 读取，但 XGL 与 WGL 尺寸相同，且只在 theta 计算中用于取节点位置
				theta.push_back(std::acos(xGL));
			}
			theta_level_k1.push_back(theta);

			for (int j = 0; j < (1 << (i + 1)) * L_k1; ++j) {
				if (j == 0) {
					WGL_phi_k1.push_back(Pi / ((1 << i) * L_k1));
				}
				phi.push_back(j * Pi / ((1 << i) * L_k1) + 0.5 * Pi / ((1 << i) * L_k1));
			}
			phi_level_k1.push_back(phi);

			std::vector<kp_Point> current_lvl_kp;
			for (int j = 0; j < static_cast<int>(theta.size()); ++j) {
				for (int k = 0; k < static_cast<int>(phi.size()); ++k) {
					double kx = std::sin(theta_level_k1[i][j]) * std::cos(phi_level_k1[i][k]);
					double ky = std::sin(theta_level_k1[i][j]) * std::sin(phi_level_k1[i][k]);
					double kz = std::cos(theta_level_k1[i][j]);
					current_lvl_kp.emplace_back(kx, ky, kz);
				}
			}
			kp_lvl_k1.push_back(current_lvl_kp);
		}

		if (integralEquType_ == 2) {
			L_k2 = static_cast<int>(std::ceil(wave.k2_abs() * D_ + 6.0 * std::cbrt(wave.k2_abs() * D_)));

			for (int i = 0; i < levelSpan; i++) {
				int n = (1 << i) * L_k2;
				AlgoFMM::gauss_legendre(n, XGL_temp, WGL_temp);
				XGL_local.push_back(XGL_temp);
				WGL_k2.push_back(WGL_temp);
				XGL_temp.clear();
				WGL_temp.clear();
			}

			for (int i = 0; i < levelSpan; i++) {
				theta.clear();
				phi.clear();

				for (int j = static_cast<int>(WGL_k2[i].size()) - 1; j >= 0; j--) {
					double xGL = XGL_local[i + levelSpan][j];  // 原来从 XGL 读取，但 XGL 与 WGL 尺寸相同，且只在 theta 计算中用于取节点位置
					theta.push_back(std::acos(xGL));
				}
				theta_level_k2.push_back(theta);

				for (int j = 0; j < (1 << (i + 1)) * L_k2; ++j) {
					if (j == 0) {
						WGL_phi_k2.push_back(Pi / ((1 << i) * L_k2));
					}
					phi.push_back(j * Pi / ((1 << i) * L_k2) + 0.5 * Pi / ((1 << i) * L_k2));
				}
				phi_level_k2.push_back(phi);

				std::vector<kp_Point> current_lvl_kp;
				for (int j = 0; j < static_cast<int>(theta.size()); ++j) {
					for (int k = 0; k < static_cast<int>(phi.size()); ++k) {
						double kx = std::sin(theta_level_k2[i][j]) * std::cos(phi_level_k2[i][k]);
						double ky = std::sin(theta_level_k2[i][j]) * std::sin(phi_level_k2[i][k]);
						double kz = std::cos(theta_level_k2[i][j]);
						current_lvl_kp.emplace_back(kx, ky, kz);
					}
				}
				kp_lvl_k2.push_back(current_lvl_kp);
			}
		}
	}
}// namespace AlgoFMM

#endif // !ALGOFMM_H
