// Zmartix.h
#pragma once
#ifndef Zmartix_H
#define Zmartix_H

#include "EFIE.h"
#include "MFIE.h"
#include "OCTree.h"
#include "EMSource.h"
#include "Operator.h"

#include <complex>
#include <vector>
#include <atomic>
#include <thread>
#include <iostream>
#include <iomanip>
#include <chrono>
#include <omp.h>

namespace Zmartix {
	inline void showProgress(int finished, int total) {
		float progress = float(finished) / total;
		std::cout << "\r[";
		int pos = int(50 * progress);
		for (int i = 0; i < 50; ++i) {
			if (i < pos) std::cout << "=";
			else if (i == pos) std::cout << ">";
			else std::cout << " ";
		}
		std::cout << "] " << std::setw(3) << int(progress * 100.0) << "%";
		std::cout.flush();
	}

	inline void computeZ_pec_efie(const std::vector<RWGBase>& rwgs, const int row,
		const EMSource& wave, const gaussPoints& gausspoint,
		std::vector<std::vector<std::complex<double>>>& Z_mom,
		bool showProgressBar = true)
	{
		auto t_start = std::chrono::high_resolution_clock::now();

		// Initialize Z_mom
		Z_mom.resize(row);
		for (int m = 0; m < row; ++m) {
			Z_mom[m].assign(row, std::complex<double>(0.0, 0.0));
		}
		std::cout << "Info: Z_mom (EFIE) initialized, starting computation... (Zmartix.h)" << std::endl;

		std::atomic<int> finishedRWG(0);
		const int totalRWG = row;
		std::thread progressThread;

		if (showProgressBar) {
			progressThread = std::thread([&finishedRWG, totalRWG]() {
				while (true) {
					int progressCount = finishedRWG.load(std::memory_order_relaxed);
					showProgress(progressCount, totalRWG);
					if (progressCount >= totalRWG) break;
					std::this_thread::sleep_for(std::chrono::milliseconds(200));
				}
				});
		}

#pragma omp parallel for schedule(dynamic)
		for (int i = 0; i < row; ++i) {
			const RWGBase* rwgField = &rwgs[i];
			for (int j = 0; j < row; ++j) {
				const RWGBase* rwgSource = &rwgs[j];
				std::complex<double> Zmn_E(0.0, 0.0);
				EFIE::computeEFIE_Zij(rwgField, rwgSource, wave.k1(), wave.eta1(), gausspoint, Zmn_E);
				Z_mom[i][j] = Zmn_E;
			}
			if (showProgressBar) {
				finishedRWG.fetch_add(1, std::memory_order_relaxed);
			}
		}
		if (showProgressBar) {
			progressThread.join();
			auto t_end = std::chrono::high_resolution_clock::now();
			auto t_elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(t_end - t_start);
			std::cout << "\nInfo: Z_mom (EFIE) done, elapsed: " << t_elapsed.count() / 1000.0
				<< " ms, total RWGs: " << totalRWG << " (Zmartix.h)" << std::endl;
		}
	}

	inline void computeZ_pec_mfie(const std::vector<RWGBase>& rwgs, const int row,
		const EMSource& wave, const gaussPoints& gausspoint,
		std::vector<std::vector<std::complex<double>>>& Z_mom,
		bool showProgressBar = true)
	{
		auto t_start = std::chrono::high_resolution_clock::now();

		// Initialize Z_mom
		Z_mom.resize(row);
		for (int m = 0; m < row; ++m) {
			Z_mom[m].assign(row, std::complex<double>(0.0, 0.0));
		}
		std::cout << "Info: Z_mom (MFIE) initialized, starting computation... (Zmartix.h)" << std::endl;

		std::atomic<int> finishedRWG(0);
		const int totalRWG = row;
		std::thread progressThread;

		if (showProgressBar) {
			progressThread = std::thread([&finishedRWG, totalRWG]() {
				while (true) {
					int progressCount = finishedRWG.load(std::memory_order_relaxed);
					showProgress(progressCount, totalRWG);
					if (progressCount >= totalRWG) break;
					std::this_thread::sleep_for(std::chrono::milliseconds(200));
				}
				});
		}

#pragma omp parallel for schedule(dynamic)
		for (int i = 0; i < row; ++i) {
			const RWGBase* rwgField = &rwgs[i];
			for (int j = 0; j < row; ++j) {
				const RWGBase* rwgSource = &rwgs[j];
				std::complex<double> Zmn_M(0.0, 0.0);
				MFIE::computeMFIE_Zij(rwgField, rwgSource, wave.k1(), gausspoint, Zmn_M);
				Z_mom[i][j] = Zmn_M;
			}
			if (showProgressBar) {
				finishedRWG.fetch_add(1, std::memory_order_relaxed);
			}
		}
		if (showProgressBar) {
			progressThread.join();
			auto t_end = std::chrono::high_resolution_clock::now();
			auto t_elapsed = std::chrono::duration_cast<std::chrono::microseconds>(t_end - t_start);
			std::cout << "\nInfo: Z_mom (MFIE) done, elapsed: " << t_elapsed.count() / 1000.0
				<< " ms, total RWGs: " << totalRWG << " (Zmartix.h)" << std::endl;
		}
	}

	inline void computeZ_pec_cfie(const std::vector<RWGBase>& rwgs, const int row,
		const EMSource& wave, const gaussPoints& gausspoint, const double alpha,
		std::vector<std::vector<std::complex<double>>>& Z_mom,
		bool showProgressBar = true)
	{
		auto t_start = std::chrono::high_resolution_clock::now();

		Z_mom.resize(row);
		for (int m = 0; m < row; ++m) {
			Z_mom[m].assign(row, std::complex<double>(0.0, 0.0));
		}
		std::cout << "Info: Z_mom (CFIE) initialized, starting computation... (Zmartix.h)" << std::endl;

		std::atomic<int> finishedRWG(0);
		const int totalRWG = row;
		std::thread progressThread;

		if (showProgressBar) {
			progressThread = std::thread([&finishedRWG, totalRWG]() {
				while (true) {
					int progressCount = finishedRWG.load(std::memory_order_relaxed);
					showProgress(progressCount, totalRWG);
					if (progressCount >= totalRWG) break;
					std::this_thread::sleep_for(std::chrono::milliseconds(200));
				}
				});
		}

#pragma omp parallel for schedule(dynamic)
		for (int i = 0; i < row; ++i) {
			const RWGBase* rwgField = &rwgs[i];
			for (int j = 0; j < row; ++j) {
				const RWGBase* rwgSource = &rwgs[j];
				std::complex<double> Zmn_E(0.0, 0.0), Zmn_M(0.0, 0.0);
				EFIE::computeEFIE_Zij(rwgField, rwgSource, wave.k1(), wave.eta1(), gausspoint, Zmn_E);
				MFIE::computeMFIE_Zij(rwgField, rwgSource, wave.k1(), gausspoint, Zmn_M);
				Z_mom[i][j] = alpha * Zmn_E + (1.0 - alpha) * wave.eta1() * Zmn_M;
			}
			if (showProgressBar) {
				finishedRWG.fetch_add(1, std::memory_order_relaxed);
			}
		}
		if (showProgressBar) {
			progressThread.join();
			auto t_end = std::chrono::high_resolution_clock::now();
			auto t_elapsed = std::chrono::duration_cast<std::chrono::microseconds>(t_end - t_start);
			std::cout << "\nInfo: Z_mom (CFIE) done, elapsed: " << t_elapsed.count() / 1000.0
				<< " ms, total RWGs: " << totalRWG << " (Zmartix.h)" << std::endl;
		}
	}

	inline void computeZ_fmm_pec_efie(
		const std::vector<OCTree::Node*>& nodes,
		const std::vector<RWGBase>& rwgs,
		std::vector<std::vector<std::complex<double>>>& Z_near,
		std::vector<std::vector<int>>& Z_near_id,
		const int row, const EMSource& wave,
		const gaussPoints& gausspoint, bool showProgressBar = true)
	{
		auto t_start = std::chrono::high_resolution_clock::now();

		Z_near.resize(row);
		Z_near_id.resize(row);
		for (int nodeIdx = 0; nodeIdx < nodes.size(); ++nodeIdx) {
			OCTree::Node* node = nodes[nodeIdx];
			int nodeRwg = node->rwgIndices.size();
			int count = nodeRwg;
			for (int nodeNR = 0; nodeNR < node->nearNeighbors.size(); ++nodeNR) {
				count += node->nearNeighbors[nodeNR]->rwgIndices.size();
			}
			for (int nodeRwgIdx = 0; nodeRwgIdx < node->rwgIndices.size(); ++nodeRwgIdx) {
				int id = node->rwgIndices[nodeRwgIdx]->rwgid;
				Z_near[id].assign(count, std::complex<double>(0.0));
				Z_near_id[id].assign(count, 0);
			}
		}

		std::atomic<int> finishedNodes(0);
		const int totalRWG = row;
		std::thread progressThread;

		if (showProgressBar) {
			progressThread = std::thread([&finishedNodes, totalRWG]() {
				while (true) {
					int progressCount = finishedNodes.load(std::memory_order_relaxed);
					showProgress(progressCount, totalRWG);
					if (progressCount >= totalRWG) break;
					std::this_thread::sleep_for(std::chrono::milliseconds(200));
				}
				});
		}

#pragma omp parallel for schedule(dynamic)
		for (int i = 0; i < nodes.size(); ++i) {
			OCTree::Node* nodeF = nodes[i];

			for (int j = 0; j < nodeF->rwgIndices.size(); ++j) {
				const int rwgfid = nodeF->rwgIndices[j]->rwgid;
				const RWGBase* rwgField = &rwgs[rwgfid];
				int count = 0;

				for (int k = 0; k < (nodeF->nearNeighbors.size() + 1); ++k) {
					if (k == nodeF->nearNeighbors.size()) {
						const OCTree::Node* nodeS = nodeF;

						for (int l = 0; l < nodeS->rwgIndices.size(); ++l) {
							int rwgSid = nodeS->rwgIndices[l]->rwgid;
							const RWGBase* rwgSource = &rwgs[rwgSid];
							std::complex<double> Zmn_near_E(0.0, 0.0);
							EFIE::computeEFIE_Zij(rwgField, rwgSource, wave.k1(), wave.eta1(), gausspoint, Zmn_near_E);
							Z_near[rwgfid][count] = Zmn_near_E;
							Z_near_id[rwgfid][count] = rwgSid;
							count++;
						}
					}
					else {
						const OCTree::Node* nodeS = nodeF->nearNeighbors[k];

						for (int l = 0; l < nodeS->rwgIndices.size(); ++l) {
							int rwgSid = nodeS->rwgIndices[l]->rwgid;
							const RWGBase* rwgSource = &rwgs[rwgSid];
							std::complex<double> Zmn_near_E(0.0, 0.0);
							EFIE::computeEFIE_Zij(rwgField, rwgSource, wave.k1(), wave.eta1(), gausspoint, Zmn_near_E);
							Z_near[rwgfid][count] = Zmn_near_E;
							Z_near_id[rwgfid][count] = rwgSid;
							count++;
						}
					}
				}
				if (showProgressBar) {
					finishedNodes.fetch_add(1, std::memory_order_relaxed);
				}
			}
		}
		if (showProgressBar) {
			progressThread.join();
			auto t_end = std::chrono::high_resolution_clock::now();
			auto t_elapsed = std::chrono::duration_cast<std::chrono::microseconds>(t_end - t_start);
			std::cout << "\nInfo: Z_near (EFIE) done, elapsed: " << t_elapsed.count() / 1000.0
				<< " ms, total RWGs: " << totalRWG << " (Zmartix.h)" << std::endl;
		}
	}

	inline void computeZ_fmm_pec_mfie(
		const std::vector<OCTree::Node*>& nodes,
		const std::vector<RWGBase>& rwgs,
		std::vector<std::vector<std::complex<double>>>& Z_near,
		std::vector<std::vector<int>>& Z_near_id,
		const int row, const EMSource& wave, const gaussPoints& gausspoint,
		bool showProgressBar = true)
	{
		auto t_start = std::chrono::high_resolution_clock::now();

		Z_near.resize(row);
		Z_near_id.resize(row);
		for (int nodeIdx = 0; nodeIdx < nodes.size(); ++nodeIdx) {
			OCTree::Node* node = nodes[nodeIdx];
			int nodeRwg = node->rwgIndices.size();
			int count = nodeRwg;
			for (int nodeNR = 0; nodeNR < node->nearNeighbors.size(); ++nodeNR) {
				count += node->nearNeighbors[nodeNR]->rwgIndices.size();
			}
			for (int nodeRwgIdx = 0; nodeRwgIdx < node->rwgIndices.size(); ++nodeRwgIdx) {
				int id = node->rwgIndices[nodeRwgIdx]->rwgid;
				Z_near[id].assign(count, std::complex<double>(0.0));
				Z_near_id[id].assign(count, 0);
			}
		}

		std::atomic<int> finishedNodes(0);
		const int totalRWG = row;
		std::thread progressThread;

		if (showProgressBar) {
			progressThread = std::thread([&finishedNodes, totalRWG]() {
				while (true) {
					int progressCount = finishedNodes.load(std::memory_order_relaxed);
					showProgress(progressCount, totalRWG);
					if (progressCount >= totalRWG) break;
					std::this_thread::sleep_for(std::chrono::milliseconds(200));
				}
				});
		}

#pragma omp parallel for schedule(dynamic)
		for (int i = 0; i < nodes.size(); ++i) {
			OCTree::Node* nodeF = nodes[i];

			for (int j = 0; j < nodeF->rwgIndices.size(); ++j) {
				const int rwgfid = nodeF->rwgIndices[j]->rwgid;
				const RWGBase* rwgField = &rwgs[rwgfid];
				int count = 0;

				for (int k = 0; k < (nodeF->nearNeighbors.size() + 1); ++k) {
					if (k == nodeF->nearNeighbors.size()) {
						const OCTree::Node* nodeS = nodeF;

						for (int l = 0; l < nodeS->rwgIndices.size(); ++l) {
							int rwgSid = nodeS->rwgIndices[l]->rwgid;
							const RWGBase* rwgSource = &rwgs[rwgSid];
							std::complex<double> Zmn_near_M(0.0, 0.0);
							MFIE::computeMFIE_Zij(rwgField, rwgSource, wave.k1(), gausspoint, Zmn_near_M);
							Z_near[rwgfid][count] = Zmn_near_M;
							Z_near_id[rwgfid][count] = rwgSid;
							count++;
						}
					}
					else {
						const OCTree::Node* nodeS = nodeF->nearNeighbors[k];

						for (int l = 0; l < nodeS->rwgIndices.size(); ++l) {
							int rwgSid = nodeS->rwgIndices[l]->rwgid;
							const RWGBase* rwgSource = &rwgs[rwgSid];
							std::complex<double> Zmn_near_M(0.0, 0.0);
							MFIE::computeMFIE_Zij(rwgField, rwgSource, wave.k1(), gausspoint, Zmn_near_M);
							Z_near[rwgfid][count] = Zmn_near_M;
							Z_near_id[rwgfid][count] = rwgSid;
							count++;
						}
					}
				}
				if (showProgressBar) {
					finishedNodes.fetch_add(1, std::memory_order_relaxed);
				}
			}
		}
		if (showProgressBar) {
			progressThread.join();
			auto t_end = std::chrono::high_resolution_clock::now();
			auto t_elapsed = std::chrono::duration_cast<std::chrono::microseconds>(t_end - t_start);
			std::cout << "\nInfo: Z_near (EFIE) done, elapsed: " << t_elapsed.count() / 1000.0
				<< " ms, total RWGs: " << totalRWG << " (Zmartix.h)" << std::endl;
		}
	}

	inline void computeZ_fmm_pec_cfie(
		const std::vector<OCTree::Node*>& nodes,
		const std::vector<RWGBase>& rwgs,
		std::vector<std::vector<std::complex<double>>>& Z_near,
		std::vector<std::vector<int>>& Z_near_id,
		const int row, const EMSource& wave, const gaussPoints& gausspoint,
		const double alpha, bool showProgressBar = true)
	{
		auto t_start = std::chrono::high_resolution_clock::now();

		Z_near.resize(row);
		Z_near_id.resize(row);
		for (int nodeIdx = 0; nodeIdx < nodes.size(); ++nodeIdx) {
			OCTree::Node* node = nodes[nodeIdx];
			int nodeRwg = node->rwgIndices.size();
			int count = nodeRwg;
			for (int nodeNR = 0; nodeNR < node->nearNeighbors.size(); ++nodeNR) {
				count += node->nearNeighbors[nodeNR]->rwgIndices.size();
			}
			for (int nodeRwgIdx = 0; nodeRwgIdx < node->rwgIndices.size(); ++nodeRwgIdx) {
				int id = node->rwgIndices[nodeRwgIdx]->rwgid;
				Z_near[id].assign(count, std::complex<double>(0.0));
				Z_near_id[id].assign(count, 0);
			}
		}

		std::atomic<int> finishedNodes(0);
		const int totalRWG = row;
		std::thread progressThread;

		if (showProgressBar) {
			progressThread = std::thread([&finishedNodes, totalRWG]() {
				while (true) {
					int progressCount = finishedNodes.load(std::memory_order_relaxed);
					showProgress(progressCount, totalRWG);
					if (progressCount >= totalRWG) break;
					std::this_thread::sleep_for(std::chrono::milliseconds(200));
				}
				});
		}
#pragma omp parallel for schedule(dynamic)
		for (int i = 0; i < nodes.size(); ++i) {
			OCTree::Node* nodeF = nodes[i];

			for (int j = 0; j < nodeF->rwgIndices.size(); ++j) {
				const int rwgfid = nodeF->rwgIndices[j]->rwgid;
				const RWGBase* rwgField = &rwgs[rwgfid];
				int count = 0;

				for (int k = 0; k < (nodeF->nearNeighbors.size() + 1); ++k) {
					if (k == nodeF->nearNeighbors.size()) {
						const OCTree::Node* nodeS = nodeF;

						for (int l = 0; l < nodeS->rwgIndices.size(); ++l) {
							int rwgSid = nodeS->rwgIndices[l]->rwgid;
							const RWGBase* rwgSource = &rwgs[rwgSid];

							std::complex<double> Zmn_near_M(0.0, 0.0), Zmn_near_E(0.0, 0.0);
							EFIE::computeEFIE_Zij(rwgField, rwgSource, wave.k1(), wave.eta1(), gausspoint, Zmn_near_E);
							MFIE::computeMFIE_Zij(rwgField, rwgSource, wave.k1(), gausspoint, Zmn_near_M);
							Z_near[rwgfid][count] = alpha * Zmn_near_E + (1.0 - alpha) * wave.eta1() * Zmn_near_M;
							Z_near_id[rwgfid][count] = rwgSid;
							count++;
						}
					}
					else {
						const OCTree::Node* nodeS = nodeF->nearNeighbors[k];

						for (int l = 0; l < nodeS->rwgIndices.size(); ++l) {
							int rwgSid = nodeS->rwgIndices[l]->rwgid;
							const RWGBase* rwgSource = &rwgs[rwgSid];

							std::complex<double> Zmn_near_M(0.0, 0.0), Zmn_near_E(0.0, 0.0);
							EFIE::computeEFIE_Zij(rwgField, rwgSource, wave.k1(), wave.eta1(), gausspoint, Zmn_near_E);
							MFIE::computeMFIE_Zij(rwgField, rwgSource, wave.k1(), gausspoint, Zmn_near_M);
							Z_near[rwgfid][count] = alpha * Zmn_near_E + (1.0 - alpha) * wave.eta1() * Zmn_near_M;
							Z_near_id[rwgfid][count] = rwgSid;
							count++;
						}
					}
				}
				if (showProgressBar) {
					finishedNodes.fetch_add(1, std::memory_order_relaxed);
				}
			}
		}
		if (showProgressBar) {
			progressThread.join();
			auto t_end = std::chrono::high_resolution_clock::now();
			auto t_elapsed = std::chrono::duration_cast<std::chrono::microseconds>(t_end - t_start);
			std::cout << "\nInfo: Z_near (CFIE) done, elapsed: " << t_elapsed.count() / 1000.0
				<< " ms, total RWGs: " << totalRWG << " (Zmartix.h)" << std::endl;
		}
	}

	inline void computeZ_die_pmchwt(
		const std::vector<RWGBase>& rwgs,
		const int row, const EMSource& wave, const gaussPoints& gausspoint,
		std::vector<std::vector<std::complex<double>>>& Z_mom,
		bool showProgressBar = true)
	{
		auto t_start = std::chrono::high_resolution_clock::now();

		const int N = row;
		const int N2 = 2 * N;
		Z_mom.resize(N2);
		for (int m = 0; m < N2; ++m) {
			Z_mom[m].assign(N2, std::complex<double>(0.0, 0.0));
		}

		std::cout << "Info: Z_mom (PMCHWT, 2N x 2N) initialized. "
			<< "N = " << N << ", k1 = " << wave.k1() << ", k2 = " << wave.k2()
			<< ", eta1=" << wave.eta1() << ", eta2=" << wave.eta2()
			<< ", epsR=" << wave.epsilonR() << ", muR=" << wave.muR() << std::endl;

		std::atomic<int> finishedRWG(0);
		const int totalRWG = N;
		std::thread progressThread;

		if (showProgressBar) {
			progressThread = std::thread([&finishedRWG, totalRWG]() {
				while (true) {
					int p = finishedRWG.load(std::memory_order_relaxed);
					showProgress(p, totalRWG);
					if (p >= totalRWG) break;
					std::this_thread::sleep_for(std::chrono::milliseconds(200));
				}
				});
		}

#pragma omp parallel for schedule(dynamic)
		for (int i = 0; i < N; ++i) {
			const RWGBase* rwgField = &rwgs[i];
			for (int j = 0; j < N; ++j) {
				const RWGBase* rwgSource = &rwgs[j];

				std::complex<double> ZL1(0.0, 0.0), ZL2(0.0, 0.0);
				std::complex<double> ZK1(0.0, 0.0), ZK2(0.0, 0.0);

				OLK::L_operator(rwgField, rwgSource, wave.k1(), wave.eta1(), gausspoint, ZL1);
				OLK::L_operator(rwgField, rwgSource, wave.k2(), wave.eta2(), gausspoint, ZL2);

				OLK::K_operator(rwgField, rwgSource, wave.k1(), gausspoint, ZK1);
				OLK::K_operator(rwgField, rwgSource, wave.k2(), gausspoint, ZK2);

				std::complex<double> z11 = ZL1 + wave.etai() * ZL2;
				std::complex<double> z12 = -(ZK1 + ZK2);
				std::complex<double> z21 = ZK1 + ZK2;
				std::complex<double> z22 = ZL1 + (1.0 / wave.etai()) * ZL2;

				Z_mom[i][j] = z11;
				Z_mom[i][N + j] = z12;
				Z_mom[N + i][j] = z21;
				Z_mom[N + i][N + j] = z22;
			}
			if (showProgressBar) {
				finishedRWG.fetch_add(1, std::memory_order_relaxed);
			}
		}
		if (showProgressBar) {
			progressThread.join();
			auto t_end = std::chrono::high_resolution_clock::now();
			auto t_elapsed = std::chrono::duration_cast<std::chrono::microseconds>(t_end - t_start);
			std::cout << "\nInfo: Z_mom (PMCHWT) done, elapsed: " << t_elapsed.count() / 1000.0
				<< " ms, total RWGs: " << totalRWG
				<< ", matrix: " << N2 << "x" << N2 << " (Zmartix.h)" << std::endl;
		}
	}

	inline void computeZ_fmm_die_pmchwt(
		const std::vector<OCTree::Node*>& nodes,
		const std::vector<RWGBase>& rwgs,
		std::vector<std::vector<std::complex<double>>>& Z_near,
		std::vector<std::vector<int>>& Z_near_id, const int row,
		const EMSource& wave, const gaussPoints& gausspoint, bool showProgressBar = true)
	{
		auto t_start = std::chrono::high_resolution_clock::now();

		const int N2 = 2 * row;
		Z_near.resize(N2);
		Z_near_id.resize(N2);

		for (int nodeIdx = 0; nodeIdx < nodes.size(); ++nodeIdx) {
			OCTree::Node* node = nodes[nodeIdx];
			int nodeRwg = node->rwgIndices.size();
			int nearCount = nodeRwg;
			for (int nr = 0; nr < node->nearNeighbors.size(); ++nr) {
				nearCount += node->nearNeighbors[nr]->rwgIndices.size();
			}

			for (int ri = 0; ri < node->rwgIndices.size(); ++ri) {
				int id = node->rwgIndices[ri]->rwgid;
				Z_near[id].assign(2 * nearCount, std::complex<double>(0.0));
				Z_near_id[id].assign(2 * nearCount, 0);
				Z_near[row + id].assign(2 * nearCount, std::complex<double>(0.0));
				Z_near_id[row + id].assign(2 * nearCount, 0);
			}
		}

		std::atomic<int> finishedRWG(0);
		const int totalRWG = row;
		std::thread progressThread;

		if (showProgressBar) {
			progressThread = std::thread([&]() {
				while (true) {
					int p = finishedRWG.load(std::memory_order_relaxed);
					showProgress(p, totalRWG);
					if (p >= totalRWG) break;
					std::this_thread::sleep_for(std::chrono::milliseconds(200));
				}
				});
		}

#pragma omp parallel for schedule(dynamic)
		for (int i = 0; i < nodes.size(); ++i) {
			OCTree::Node* nodeF = nodes[i];

			for (int j = 0; j < nodeF->rwgIndices.size(); ++j) {
				const int rwgfid = nodeF->rwgIndices[j]->rwgid;
				const RWGBase* rwgField = &rwgs[rwgfid];
				int count = 0;

				for (int k = 0; k < (nodeF->nearNeighbors.size() + 1); ++k) {
					const OCTree::Node* nodeS = (k == nodeF->nearNeighbors.size()) ? nodeF : nodeF->nearNeighbors[k];

					for (int l = 0; l < nodeS->rwgIndices.size(); ++l) {
						int rwgSid = nodeS->rwgIndices[l]->rwgid;
						const RWGBase* rwgSource = &rwgs[rwgSid];

						std::complex<double> ZL1(0.0, 0.0), ZL2(0.0, 0.0);
						std::complex<double> ZK1(0.0, 0.0), ZK2(0.0, 0.0);

						EFIE::computeEFIE_Zij(rwgField, rwgSource, wave.k1(), wave.eta1(), gausspoint, ZL1);
						EFIE::computeEFIE_Zij(rwgField, rwgSource, wave.k2(), wave.eta2(), gausspoint, ZL2);

						MFIE::computeMFIE_Zij(rwgField, rwgSource, wave.k1(), gausspoint, ZK1);
						MFIE::computeMFIE_Zij(rwgField, rwgSource, wave.k2(), gausspoint, ZK2);

						Z_near[rwgfid][count] = ZL1 + ZL2;
						Z_near_id[rwgfid][count] = rwgSid;
						count++;

						Z_near[rwgfid][count] = -ZK1 - ZK2;
						Z_near_id[rwgfid][count] = row + rwgSid;
						count++;

						Z_near[row + rwgfid][count - 2] = ZK1 + ZK2;
						Z_near[row + rwgfid][count - 1] = (std::complex<double>(1, 0) / pow(wave.eta1(), 2)) * ZL1 +
							(std::complex<double>(1, 0) / pow(wave.eta2(), 2)) * ZL2;

						Z_near_id[row + rwgfid][count - 2] = rwgSid;
						Z_near_id[row + rwgfid][count - 1] = row + rwgSid;
					}
				}
				if (showProgressBar) {
					finishedRWG.fetch_add(1, std::memory_order_relaxed);
				}
			}
		}
		if (showProgressBar) {
			progressThread.join();
			auto t_end = std::chrono::high_resolution_clock::now();
			auto t_elapsed = std::chrono::duration_cast<std::chrono::microseconds>(t_end - t_start);
			std::cout << "\nInfo: Z_near (PMCHWT) done, elapsed: " << t_elapsed.count() / 1000.0
				<< " ms, total RWGs: " << totalRWG << " (Zmartix.h)" << std::endl;
		}
	}
}

#endif // !Zmartix_H
