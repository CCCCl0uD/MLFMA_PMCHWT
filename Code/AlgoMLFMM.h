#pragma once
#ifndef ALGOMLFMM_H
#define ALGOMLFMM_H

#include <atomic>
#include <algorithm>
#include <vector>
#include <chrono>
#include <complex>
#include <array>
#include <thread>
#include <iostream>
#include <iomanip>

#include "OverloadAlgo.h"
#include "OCTree.h"
#include "Bessel.h"
#include "EMSource.h"

namespace AlgoMLFMM {
	inline void computeBase_MLFMM(
		const std::vector<std::vector<OCTree::Node*>>& octreeNodes,
		const EMSource& wave, const int integralEquType_,
		int maxLevel_, int rwgSize, std::vector<int>& L_k1, std::vector<int>& L_k2, int& row, int& levelSpan,
		std::vector<std::vector<double>>& WGL_k1, std::vector<std::vector<double>>& WGL_k2,
		std::vector<double>& WGL_phi_k1, std::vector<double>& WGL_phi_k2,
		std::vector<std::vector<double>>& theta_level_k1, std::vector<std::vector<double>>& theta_level_k2,
		std::vector<std::vector<double>>& phi_level_k1, std::vector<std::vector<double>>& phi_level_k2,
		std::vector<std::vector<kp_Point>>& kp_lvl_k1, std::vector<std::vector<kp_Point>>& kp_lvl_k2)
	{
		row = rwgSize;
		levelSpan = maxLevel_ - 1;

		auto cal_L_k = [&](double kAbs, int level)->int {
			const int physicalLevel = maxLevel_ - level; // 物理层级
			const double highLevelLen = octreeNodes[physicalLevel][0]->cubeLength;
			const double D_ = std::sqrt(3.0) * highLevelLen;
			return static_cast<int>(std::ceil(kAbs * D_ + 6.0 * std::cbrt(kAbs * D_)));
			};

		auto buildAngularGrid = [&](double kAbs,
			std::vector<int>& L_by_level,
			std::vector<std::vector<double>>& WGL,
			std::vector<double>& WGL_phi,
			std::vector<std::vector<double>>& theta_level,
			std::vector<std::vector<double>>& phi_level,
			std::vector<std::vector<kp_Point>>& kp_lvl) {
				L_by_level.clear();
				WGL.clear();
				WGL_phi.clear();
				theta_level.clear();
				phi_level.clear();
				kp_lvl.clear();

				L_by_level.reserve(levelSpan);
				WGL.reserve(levelSpan);
				WGL_phi.reserve(levelSpan);
				theta_level.reserve(levelSpan);
				phi_level.reserve(levelSpan);
				kp_lvl.reserve(levelSpan);

				for (int i = 0; i < levelSpan; ++i) {
					const int L = cal_L_k(kAbs, i);
					const int thetaCount = L;
					//const int thetaCount = L + 1;

					L_by_level.push_back(L);

					std::vector<double> XGL_temp;
					std::vector<double> WGL_temp;
					my_math::gauss_legendre(thetaCount, XGL_temp, WGL_temp);
					WGL.push_back(WGL_temp);

					std::vector<double> theta;
					theta.reserve(thetaCount);
					for (int j = thetaCount - 1; j >= 0; --j) {
						theta.push_back(std::acos(XGL_temp[j]));
					}
					theta_level.push_back(theta);

					std::vector<double> phi;
					const int phiCount = 2 * L;
					const double dphi = Pi / L;
					phi.reserve(phiCount);
					WGL_phi.push_back(dphi);
					for (int j = 0; j < phiCount; ++j) {
						phi.push_back(j * dphi + 0.5 * dphi);
					}
					phi_level.push_back(phi);

					std::vector<kp_Point> current_lvl_kp;
					current_lvl_kp.reserve(theta_level[i].size() * phi_level[i].size());
					for (int t = 0; t < static_cast<int>(theta_level[i].size()); ++t) {
						for (int p = 0; p < static_cast<int>(phi_level[i].size()); ++p) {
							double kx = std::sin(theta_level[i][t]) * std::cos(phi_level[i][p]);
							double ky = std::sin(theta_level[i][t]) * std::sin(phi_level[i][p]);
							double kz = std::cos(theta_level[i][t]);
							current_lvl_kp.emplace_back(kx, ky, kz);
						}
					}
					kp_lvl.push_back(current_lvl_kp);
				}
			};

		buildAngularGrid(wave.k1_abs(), L_k1, WGL_k1, WGL_phi_k1, theta_level_k1, phi_level_k1, kp_lvl_k1);

		if (integralEquType_ == 2) {
			buildAngularGrid(wave.k2_abs(), L_k2, WGL_k2, WGL_phi_k2, theta_level_k2, phi_level_k2, kp_lvl_k2);
		}

		auto printAngularGridInfo = [&](const char* tag,
			const std::vector<int>& L_by_level,
			const std::vector<std::vector<double>>& theta_level,
			const std::vector<std::vector<double>>& phi_level,
			const std::vector<std::vector<kp_Point>>& kp_lvl) {
				std::cout << "Info: computeBase_MLFMM " << tag
					<< " angular grid summary (AlgoMLFMM.h)" << std::endl;
				for (int i = 0; i < static_cast<int>(L_by_level.size()); ++i) {
					const int physicalLevel = maxLevel_ - i;
					std::cout << "  physical level " << std::setw(2) << physicalLevel
						<< ": L = " << std::setw(4) << L_by_level[i]
						<< ", theta = " << std::setw(5) << theta_level[i].size()
						<< ", phi = " << std::setw(5) << phi_level[i].size()
						<< ", kp = " << std::setw(8) << kp_lvl[i].size()
						<< std::endl;
				}
			};

		std::cout << "Info: computeBase_MLFMM done, row = " << row
			<< ", maxLevel = " << maxLevel_
			<< ", levelSpan = " << levelSpan
			<< " (AlgoMLFMM.h)" << std::endl;
		printAngularGridInfo("k1", L_k1, theta_level_k1, phi_level_k1, kp_lvl_k1);

		if (integralEquType_ == 2) {
			printAngularGridInfo("k2", L_k2, theta_level_k2, phi_level_k2, kp_lvl_k2);
		}
	}

	// 2. 插值位置查找（二分搜索）
	inline std::array<int, 3> find_interp_pos(double target, const std::vector<double>& arr, bool periodic = true) {
		int n = static_cast<int>(arr.size());
		if (n < 3) return { -1, -1, -1 };
		double delta = arr[1] - arr[0];
		bool ascending = delta > 0;
		int lo = 0, hi = n - 1, mid = 0;
		while (lo <= hi) {
			mid = (lo + hi) / 2;
			if (std::abs(arr[mid] - target) < 1e-12) break;
			if (ascending) {
				if (arr[mid] > target) hi = mid - 1;
				else lo = mid + 1;
			}
			else {
				if (arr[mid] > target) lo = mid + 1;
				else hi = mid - 1;
			}
		}
		int closest = mid;
		if (mid > 0 && std::abs(arr[mid - 1] - target) < std::abs(arr[closest] - target))
			closest = mid - 1;
		if (mid < n - 1 && std::abs(arr[mid + 1] - target) < std::abs(arr[closest] - target))
			closest = mid + 1;
		int prev, curr, next;
		curr = closest;
		if (ascending) {
			if (closest == 0) {
				if (periodic) { prev = n - 1; next = closest + 1; }
				else { prev = 0; next = 1; }
			}
			else if (closest == n - 1) {
				if (periodic) { prev = closest - 1; next = 0; }
				else { prev = closest - 1; next = closest; }
			}
			else { prev = closest - 1; next = closest + 1; }
		}
		else {
			if (closest == 0) {
				if (periodic) { prev = closest + 1; next = n - 1; }
				else { prev = 0; next = 1; }
			}
			else if (closest == n - 1) {
				if (periodic) { prev = 0; next = closest - 1; }
				else { prev = closest - 1; next = closest; }
			}
			else { prev = closest + 1; next = closest - 1; }
		}
		return { prev, curr, next };
	}

	// 3. 计算转移因子矩阵
	inline void computeTransferFactorMatrix(
		std::vector<std::vector<std::vector<std::complex<double>>>>& TF,
		const std::vector<std::vector<OCTree::nodePoint>>& octreeNodesDRvec,
		const std::vector<std::vector<kp_Point>>& kp_lvl,
		const std::vector<int>& L_lvl,
		int maxLevel_, std::complex<double> k_wavenumber)
	{
		auto t_start = std::chrono::high_resolution_clock::now();

		const int span = maxLevel_ - 1;

		if (maxLevel_ - 1 <= 0) {
			std::cerr << "Warning: maxLevel_ too small, cannot compute transfer matrix" << std::endl;
			return;
		}
		if (kp_lvl.empty()) {
			std::cerr << "Error: kp_level not initialized" << std::endl;
			return;
		}
		if (static_cast<int>(kp_lvl.size()) < span || static_cast<int>(L_lvl.size()) < span) {
			std::cerr << "Error: kp_lvl or L_lvl size mismatch in computeTransferFactorMatrix" << std::endl;
			return;
		}
		TF.resize(maxLevel_ - 1);
		for (int i = 0; i < (maxLevel_ - 1); ++i) {
			int vec_at_level = static_cast<int>(octreeNodesDRvec[i].size());
			TF[i].resize(vec_at_level);
			int rev = span - 1 - i;
			int kp_at_level = static_cast<int>(kp_lvl[rev].size());

			for (int j = 0; j < vec_at_level; ++j) {
				TF[i][j].assign(kp_at_level, std::complex<double>(0.0));
			}
		}
		for (int level = 0; level < maxLevel_ - 1; ++level) {
			int totalVec = static_cast<int>(octreeNodesDRvec[level].size());
			int rev = span - 1 - level;

			const int L_ = L_lvl[rev];
			const int kps = static_cast<int>(kp_lvl[rev].size());

			std::atomic<int> finishedVec(0);
			const int barWidth = 50;
			std::thread progressThread([&]() {
				while (true) {
					int progressCount = finishedVec.load(std::memory_order_relaxed);
					float progress = float(progressCount) / totalVec;
					std::cout << "\rLevel " << std::setw(2) << level + 2 << " [";
					int pos = int(barWidth * progress);
					for (int i = 0; i < barWidth; ++i) {
						if (i < pos) std::cout << "=";
						else if (i == pos) std::cout << ">";
						else std::cout << " ";
					}
					std::cout << "] " << std::setw(3) << int(progress * 100.0) << "%";
					std::cout.flush();
					if (progressCount >= totalVec) break;
					std::this_thread::sleep_for(std::chrono::milliseconds(300));
				}
				});
#pragma omp parallel for schedule(dynamic)
			for (int vec = 0; vec < totalVec; ++vec) {
				double r[3] = {
					octreeNodesDRvec[level][vec].x,
					octreeNodesDRvec[level][vec].y,
					octreeNodesDRvec[level][vec].z
				};
				double Rpq = std::sqrt(r[0] * r[0] + r[1] * r[1] + r[2] * r[2]);
				my_math::normalize(r, Rpq);

				std::complex<double> kR = k_wavenumber * Rpq;
				int L = std::min(L_, static_cast<int>(std::ceil(std::abs(kR))));
				std::vector<std::complex<double>> h2_array;
				Bessel::spherical_hankel_2_array(L, std::abs(kR), h2_array);

				auto& kpLevel = kp_lvl[rev];

				for (int kpNum = 0; kpNum < kps; ++kpNum) {
					double kp[3] = {
						kpLevel[kpNum].x, kpLevel[kpNum].y, kpLevel[kpNum].z
					};
					double r_dot_kp = my_math::legendrePolynomialParameters(r, kp);

					std::complex<double> sum(0.0, 0.0);
					for (int l = 0; l <= L; ++l) {
						double p = boost::math::legendre_p(l, r_dot_kp);
						auto jl = my_math::jl_table[l & 3];
						sum += jl * h2_array[l] * (2.0 * l + 1.0) * p;
					}
					TF[level][vec][kpNum] = sum;
				}
				finishedVec.fetch_add(1, std::memory_order_relaxed);
			}
			progressThread.join();
			std::cout << "\nInfo: Level " << level + 2
				<< " transfer factor done, total vectors: " << totalVec << " (AlgoMLFMM.h)" << std::endl;
			std::cout << "Info: Level " << level + 2
				<< " transfer factor done, kp points: " << kps << ", L = " << L_ << " (AlgoMLFMM.h)" << std::endl;
		}
		auto t_end = std::chrono::high_resolution_clock::now();
		auto t_elapsed = std::chrono::duration_cast<std::chrono::microseconds>(t_end - t_start);
		std::cout << "Info: Transfer factor matrix done, elapsed: "
			<< t_elapsed.count() / 1000.0 << " ms (AlgoMLFMM.h)\n" << std::endl;
	}

	// 4. 计算插值矩阵
	inline void computeInterpolationMatrix(
		std::vector<std::vector<InterpData>>& interpol,
		const std::vector<std::vector<kp_Point>>& kp_lvl,
		const std::vector<std::vector<double>>& theta_level,
		const std::vector<std::vector<double>>& phi_level,
		int maxLevel_)
	{
		auto t_start = std::chrono::high_resolution_clock::now();
		const int span = maxLevel_ - 1;
		int levels = maxLevel_ - 2;
		interpol.resize(levels);
		for (int lvl = 0; lvl < levels; ++lvl) {
			int rev = span - 1 - lvl;
			interpol[lvl].resize(kp_lvl[rev].size());
		}
		for (int level = 0; level < levels; ++level) {
			const int revF = span - 1 - level;
			const int revS = span - 2 - level;
			int tfn = static_cast<int>(theta_level[revF].size());
			int pfn = static_cast<int>(phi_level[revF].size());
			int psn = static_cast<int>(phi_level[revS].size());

#pragma omp parallel for collapse(2) schedule(static)
			for (int t = 0; t < tfn; ++t) {
				for (int p = 0; p < pfn; ++p) {
					double theta0 = theta_level[revF][t];
					auto theta_idx = find_interp_pos(theta0, theta_level[revS], true);
					double theta_vals[3] = {
						theta_level[revS][theta_idx[0]],
						theta_level[revS][theta_idx[1]],
						theta_level[revS][theta_idx[2]]
					};
					int th1_base = theta_idx[0] * psn;
					int th2_base = theta_idx[1] * psn;
					int th3_base = theta_idx[2] * psn;
					double l1T = my_math::lagrangeWeight(theta0, theta_vals[0], theta_vals[1], theta_vals[2]);
					double l2T = my_math::lagrangeWeight(theta0, theta_vals[1], theta_vals[0], theta_vals[2]);
					double l3T = my_math::lagrangeWeight(theta0, theta_vals[2], theta_vals[0], theta_vals[1]);
					double phi0 = phi_level[revF][p];
					auto phi_idx = find_interp_pos(phi0, phi_level[revS], true);
					double phi_vals[3] = {
						phi_level[revS][phi_idx[0]],
						phi_level[revS][phi_idx[1]],
						phi_level[revS][phi_idx[2]]
					};
					double l1P = my_math::lagrangeWeight(phi0, phi_vals[0], phi_vals[1], phi_vals[2]);
					double l2P = my_math::lagrangeWeight(phi0, phi_vals[1], phi_vals[0], phi_vals[2]);
					double l3P = my_math::lagrangeWeight(phi0, phi_vals[2], phi_vals[0], phi_vals[1]);
					int idx = t * pfn + p;
					std::array<double, 9> weight = {
						l1T * l1P, l1T * l2P, l1T * l3P,
						l2T * l1P, l2T * l2P, l2T * l3P,
						l3T * l1P, l3T * l2P, l3T * l3P
					};
					interpol[level][idx] = InterpData(
						idx,
						std::array<int, 3>{theta_idx[0], theta_idx[1], theta_idx[2]},
						std::array<int, 3>{phi_idx[0], phi_idx[1], phi_idx[2]},
						weight, psn);
				}
			}
		}
		for (int i = 0; i < maxLevel_ - 2; ++i) {
			std::cout << "Info: Level " << i + 2 << " interpolation size = "
				<< interpol[i].size() << " (AlgoMLFMM.h)" << std::endl;
		}
		auto t_end = std::chrono::high_resolution_clock::now();
		auto t_elapsed = std::chrono::duration_cast<std::chrono::microseconds>(t_end - t_start);
		std::cout << "Info: Interpolation matrix done, elapsed: "
			<< t_elapsed.count() / 1000.0 << " ms (AlgoMLFMM.h)\n" << std::endl;
	}
}// namespace AlgoMLFMM

#endif // !ALGOMLFMM_H
