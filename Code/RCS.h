#pragma once

#include "RCSExportConfig.h"
#include "EMSource.h"

#include <cmath>
#include <iostream>

namespace RCSUtils {
	template<typename SolverType>
	void exportI(SolverType& solver, const RCSExportConfig& cfg, const std::complex<double>* I_solve) {
		if (!I_solve) {
			std::cerr << "Error: I_solve is null, cannot export.\n";
			return;
		}
		auto iFile = cfg.openFile("_I", "txt");
		if (!iFile.is_open()) {
			std::cerr << "Error: Failed to open file for current coefficients export.\n";
			return;
		}
		iFile << "I_solve\n";
		for (int i = 0; i < solver.row; ++i) {
			iFile << i << "\t" << I_solve[i] << "\n";
		}
		iFile.close();
		std::cout << "Current coefficients export completed.\n";
	}

	template<typename SolverType>
	void exportV(SolverType& solver, const RCSExportConfig& cfg) {
		auto vFile = cfg.openFile("_V", "txt");
		if (!vFile.is_open()) {
			std::cerr << "Error: Failed to open file for voltage vector export.\n";
			return;
		}
		vFile << "RHS\n";
		for (int i = 0; i < solver.row; ++i) {
			vFile << i << "\t" << solver.Vm[i] << "\n";
		}
		vFile.close();
		std::cout << "Voltage vector export completed.\n";
	}
	inline int calculateAngleCount(double start, double end, double step) {
		const double diff = std::abs(end - start);
		if (diff < 1.0e-9) {
			return 1;
		}
		return static_cast<int>(std::ceil(diff / step - 1.0e-9)) + 1;
	}

	inline double calculateAngleAt(double start, double end, double step, int idx, int count) {
		if (idx == count - 1) {
			return end;
		}
		const double direction = (end >= start) ? 1.0 : -1.0;
		return start + direction * idx * step;
	}

	inline int calculateNumPoints(double start, double startPhi, double end, double endPhi, double step) {
		return calculateAngleCount(start, end, step) *
			calculateAngleCount(startPhi, endPhi, step);
	}

	template<typename SolverType>
	inline void calculateRCS(SolverType& solver, const std::complex<double>* I_solve, std::vector<std::complex<double>>& RCS_pre,
		std::complex<double> k_, std::complex<double> eta_, double th, double ph) {
		std::array<double, 3> r;
		const double theta = (th / 180.0) * Pi;
		const double phi = (ph / 180.0) * Pi;
		r = { sin(theta) * cos(phi), sin(theta) * sin(phi), cos(theta) };

		RCS_pre.assign(3, 0.0);

		for (int rwgid = 0; rwgid < solver.row; ++rwgid) {
			const RWGBase* rwg = &solver.rwgs[rwgid];
			double ln = rwg->edgeLength;
			std::array<std::complex<double>, 3> F_temp = { 0.0, 0.0, 0.0 };

			for (int i = 0; i < solver.gausspoint.N_points_; ++i) {
				double wf = solver.gausspoint.weight[i];

				// 正面贡献
				std::array<double, 3> pos = {
					solver.gausspoint.l1[i] * rwg->positiveFace->vertex1->x + solver.gausspoint.l2[i] * rwg->positiveFace->vertex2->x + solver.gausspoint.l3[i] * rwg->positiveFace->vertex3->x,
					solver.gausspoint.l1[i] * rwg->positiveFace->vertex1->y + solver.gausspoint.l2[i] * rwg->positiveFace->vertex2->y + solver.gausspoint.l3[i] * rwg->positiveFace->vertex3->y,
					solver.gausspoint.l1[i] * rwg->positiveFace->vertex1->z + solver.gausspoint.l2[i] * rwg->positiveFace->vertex2->z + solver.gausspoint.l3[i] * rwg->positiveFace->vertex3->z
				};
				std::array<double, 3> rho = {
					pos[0] - rwg->freeVertexPositive->x,
					pos[1] - rwg->freeVertexPositive->y,
					pos[2] - rwg->freeVertexPositive->z
				};
				std::complex<double> phase = k_ * (pos[0] * r[0] + pos[1] * r[1] + pos[2] * r[2]);
				std::complex<double> e_ = exp(J * phase);
				std::array<std::complex<double>, 3> F_pos;
				my_math::computeTensorDotVector(F_pos, r, rho);
				F_temp[0] += ln * wf * F_pos[0] * e_;
				F_temp[1] += ln * wf * F_pos[1] * e_;
				F_temp[2] += ln * wf * F_pos[2] * e_;

				// 负面贡献
				std::array<double, 3> neg = {
					solver.gausspoint.l1[i] * rwg->negativeFace->vertex1->x + solver.gausspoint.l2[i] * rwg->negativeFace->vertex2->x + solver.gausspoint.l3[i] * rwg->negativeFace->vertex3->x,
					solver.gausspoint.l1[i] * rwg->negativeFace->vertex1->y + solver.gausspoint.l2[i] * rwg->negativeFace->vertex2->y + solver.gausspoint.l3[i] * rwg->negativeFace->vertex3->y,
					solver.gausspoint.l1[i] * rwg->negativeFace->vertex1->z + solver.gausspoint.l2[i] * rwg->negativeFace->vertex2->z + solver.gausspoint.l3[i] * rwg->negativeFace->vertex3->z
				};
				rho = {
					rwg->freeVertexNegative->x - neg[0],
					rwg->freeVertexNegative->y - neg[1],
					rwg->freeVertexNegative->z - neg[2]
				};
				phase = k_ * (neg[0] * r[0] + neg[1] * r[1] + neg[2] * r[2]);
				e_ = exp(J * phase);
				std::array<std::complex<double>, 3> F_neg;
				my_math::computeTensorDotVector(F_neg, r, rho);
				F_temp[0] += ln * wf * F_neg[0] * e_;
				F_temp[1] += ln * wf * F_neg[1] * e_;
				F_temp[2] += ln * wf * F_neg[2] * e_;
			}
			std::complex<double> coeff = J * k_ * eta_;
			RCS_pre[0] += coeff * F_temp[0] * I_solve[rwgid];
			RCS_pre[1] += coeff * F_temp[1] * I_solve[rwgid];
			RCS_pre[2] += coeff * F_temp[2] * I_solve[rwgid];
		}
	}

	template<typename SolverType>
	inline void calculateRCS_PMCHWT(SolverType& solver, const std::complex<double>* I_J, const std::complex<double>* I_M,
		std::vector<std::complex<double>>& RCS_pre, const std::complex<double> k1, const std::complex<double> k2,
		const std::complex<double> eta1, const std::complex<double> eta2, double th, double ph)
	{
		const double theta = (th / 180.0) * Pi;
		const double phi = (ph / 180.0) * Pi;
		std::array<double, 3> r = { sin(theta) * cos(phi), sin(theta) * sin(phi), cos(theta) };

		RCS_pre.assign(3, 0.0);

		for (int rwgid = 0; rwgid < solver.row; ++rwgid) {
			const RWGBase* rwg = &solver.rwgs[rwgid];
			double ln = rwg->edgeLength;
			std::array<std::complex<double>, 3> F_J = { 0.0,0.0,0.0 };
			std::array<std::complex<double>, 3> F_M = { 0.0,0.0,0.0 };

			for (int i = 0; i < solver.gausspoint.N_points_; ++i) {
				double wf = solver.gausspoint.weight[i];

				std::array<double, 3> pos = {
					solver.gausspoint.l1[i] * rwg->positiveFace->vertex1->x + solver.gausspoint.l2[i] * rwg->positiveFace->vertex2->x + solver.gausspoint.l3[i] * rwg->positiveFace->vertex3->x,
					solver.gausspoint.l1[i] * rwg->positiveFace->vertex1->y + solver.gausspoint.l2[i] * rwg->positiveFace->vertex2->y + solver.gausspoint.l3[i] * rwg->positiveFace->vertex3->y,
					solver.gausspoint.l1[i] * rwg->positiveFace->vertex1->z + solver.gausspoint.l2[i] * rwg->positiveFace->vertex2->z + solver.gausspoint.l3[i] * rwg->positiveFace->vertex3->z
				};

				std::array<double, 3> rho = {
					pos[0] - rwg->freeVertexPositive->x,
					pos[1] - rwg->freeVertexPositive->y,
					pos[2] - rwg->freeVertexPositive->z
				};
				std::complex<double> phase1 = k1 * (pos[0] * r[0] + pos[1] * r[1] + pos[2] * r[2]);
				std::complex<double> e_k1 = exp(J * phase1);

				std::array<std::complex<double>, 3> Fp;
				my_math::computeTensorDotVector(Fp, r, rho);
				F_J[0] += ln * wf * Fp[0] * e_k1;
				F_J[1] += ln * wf * Fp[1] * e_k1;
				F_J[2] += ln * wf * Fp[2] * e_k1;

				std::array<double, 3> rCrossRho = cross(r, rho);
				F_M[0] += ln * wf * rCrossRho[0] * e_k1;
				F_M[1] += ln * wf * rCrossRho[1] * e_k1;
				F_M[2] += ln * wf * rCrossRho[2] * e_k1;

				std::array<double, 3> neg = {
					solver.gausspoint.l1[i] * rwg->negativeFace->vertex1->x + solver.gausspoint.l2[i] * rwg->negativeFace->vertex2->x + solver.gausspoint.l3[i] * rwg->negativeFace->vertex3->x,
					solver.gausspoint.l1[i] * rwg->negativeFace->vertex1->y + solver.gausspoint.l2[i] * rwg->negativeFace->vertex2->y + solver.gausspoint.l3[i] * rwg->negativeFace->vertex3->y,
					solver.gausspoint.l1[i] * rwg->negativeFace->vertex1->z + solver.gausspoint.l2[i] * rwg->negativeFace->vertex2->z + solver.gausspoint.l3[i] * rwg->negativeFace->vertex3->z
				};

				rho = {
					rwg->freeVertexNegative->x - neg[0],
					rwg->freeVertexNegative->y - neg[1],
					rwg->freeVertexNegative->z - neg[2]
				};
				phase1 = k1 * (neg[0] * r[0] + neg[1] * r[1] + neg[2] * r[2]);
				e_k1 = exp(J * phase1);

				my_math::computeTensorDotVector(Fp, r, rho);
				F_J[0] += ln * wf * Fp[0] * e_k1;
				F_J[1] += ln * wf * Fp[1] * e_k1;
				F_J[2] += ln * wf * Fp[2] * e_k1;

				rCrossRho = cross(r, rho);
				F_M[0] += ln * wf * rCrossRho[0] * e_k1;
				F_M[1] += ln * wf * rCrossRho[1] * e_k1;
				F_M[2] += ln * wf * rCrossRho[2] * e_k1;
			}

			std::complex<double> coeff = -J * k1;

			RCS_pre[0] += -coeff * F_J[0] * I_J[rwgid] + coeff * F_M[0] * I_M[rwgid];
			RCS_pre[1] += -coeff * F_J[1] * I_J[rwgid] + coeff * F_M[1] * I_M[rwgid];
			RCS_pre[2] += -coeff * F_J[2] * I_J[rwgid] + coeff * F_M[2] * I_M[rwgid];
		}
	}

	template<typename SolverType, typename ComputeVFunc>
	inline void computeDualStatic(SolverType& solver, const RCSExportConfig& cfg, ComputeVFunc computeV) {
		// 从 cfg 获取参数
		double inc_th = cfg.incTheta;
		double inc_ph = cfg.incPhi;
		std::string pol_wave = cfg.polWaveStr(); // "V" or "H"
		double data_pol = (pol_wave == "V") ? 0.0 : 90.0;

		// 构造入射波
		double kinc[3], einc[3], hinc[3];
		EMSource pw;
		pw.initPW(data_pol, inc_th, inc_ph);
		kinc[0] = pw.vW(0); kinc[1] = pw.vW(1); kinc[2] = pw.vW(2);
		einc[0] = pw.vE(0); einc[1] = pw.vE(1); einc[2] = pw.vE(2);
		hinc[0] = pw.vH(0); hinc[1] = pw.vH(1); hinc[2] = pw.vH(2);

		// 计算右边向量
		computeV(solver, kinc, einc, hinc);

		// 导出 V
		//exportV(solver, cfg);

		// 求解
		int itr_max = solver.row, mr = 50;
		double tol_abs = 1.0e-8;
		double tol_rel = 1.0e-4;
		std::complex<double>* I_solve = new std::complex<double>[solver.row];
		std::complex<double>* Vec_R = new std::complex<double>[solver.row];
		for (int i = 0; i < solver.row; ++i) {
			I_solve[i] = std::complex<double>(0.0, 0.0);
			Vec_R[i] = solver.Vm[i];
		}
		solver.matrix_solver(solver.row, I_solve, Vec_R, itr_max, mr, tol_abs, tol_rel);

		// 导出 I
		//exportI(solver, cfg, I_solve);

		// 计算双站 RCS
		auto rcsFile = cfg.openFile("_RCS", "txt");
		if (!rcsFile.is_open()) throw std::runtime_error("Cannot open RCS file");
		const int thetaNum = calculateAngleCount(cfg.scaThetaStart, cfg.scaThetaEnd, cfg.scaStep);
		const int phiNum = calculateAngleCount(cfg.scaPhiStart, cfg.scaPhiEnd, cfg.scaStep);
		const int totalPoints = thetaNum * phiNum;
		int idx = 0;
		for (int it = 0; it < thetaNum; ++it) {
			double th_dual = calculateAngleAt(
				cfg.scaThetaStart, cfg.scaThetaEnd, cfg.scaStep, it, thetaNum);
			for (int ip = 0; ip < phiNum; ++ip, ++idx) {
				double ph_dual = calculateAngleAt(
					cfg.scaPhiStart, cfg.scaPhiEnd, cfg.scaStep, ip, phiNum);
				std::vector<std::complex<double>> RCS_pre(3, 0.0);
				RCSUtils::calculateRCS(solver, I_solve, RCS_pre, solver.wave.k1(), solver.wave.eta1(), th_dual, ph_dual);

				// 计算极化分量
				pw.initPW(data_pol, th_dual, ph_dual);
				double vec_v[3] = { pw.vTh(0), pw.vTh(1), pw.vTh(2) };
				double vec_h[3] = { pw.vPh(0), pw.vPh(1), pw.vPh(2) };
				std::complex<double> E_v = vec_v[0] * RCS_pre[0] + vec_v[1] * RCS_pre[1] + vec_v[2] * RCS_pre[2];
				std::complex<double> E_h = vec_h[0] * RCS_pre[0] + vec_h[1] * RCS_pre[1] + vec_h[2] * RCS_pre[2];
				double rcs = 10.0 * log10((norm(E_v) + norm(E_h)) / (4.0 * Pi));
				double rcs_v = 10.0 * log10(norm(E_v) / (4.0 * Pi));
				double rcs_h = 10.0 * log10(norm(E_h) / (4.0 * Pi));
				rcsFile << std::fixed << std::setprecision(12);
				rcsFile << std::setw(16) << th_dual << "\t" << ph_dual << "\t" << idx << "\t" << rcs << "\t" << rcs_v << "\t" << rcs_h << std::endl;
				std::cout << std::setw(16) << th_dual << "\t" << ph_dual << "\t" << idx << "\t" << rcs << "\t" << rcs_v << "\t" << rcs_h << std::endl;
			}
		}
		size_t mem = solver.computeMem();
		delete[] I_solve;
		delete[] Vec_R;
		rcsFile << "\n\nAlgorithm memory = " << mem / (1024.0 * 1024.0) << " MB\n";
		rcsFile.close();
		std::cout << "Algorithm memory = " << mem / (1024.0 * 1024.0) << " MB" << std::endl;
		std::cout << "RCS export completed." << std::endl;
	}

	template<typename SolverType, typename ComputeVFunc>
	inline void computeMonoStatic(SolverType& solver, const RCSExportConfig& cfg, ComputeVFunc computeV) {
		auto rcsFile = cfg.openFile("_RCS", "txt");
		if (!rcsFile.is_open()) throw std::runtime_error("Cannot open RCS file");
		//auto iFile = cfg.openFile("_I", "txt");
		//auto vFile = cfg.openFile("_V", "txt");
		const int thetaNum = calculateAngleCount(cfg.scaThetaStart, cfg.scaThetaEnd, cfg.scaStep);
		const int phiNum = calculateAngleCount(cfg.scaPhiStart, cfg.scaPhiEnd, cfg.scaStep);
		const int totalPoints = thetaNum * phiNum;
		int idx = 0;
		for (int it = 0; it < thetaNum; ++it) {
			double th_mono = calculateAngleAt(
				cfg.scaThetaStart, cfg.scaThetaEnd, cfg.scaStep, it, thetaNum);
			for (int ip = 0; ip < phiNum; ++ip, ++idx) {
				double ph_mono = calculateAngleAt(
					cfg.scaPhiStart, cfg.scaPhiEnd, cfg.scaStep, ip, phiNum);
				std::string pol_wave = cfg.polWaveStr();
				double data_pol = (pol_wave == "V") ? 0.0 : 90.0;

				double kinc[3], einc[3], hinc[3];
				EMSource pw;
				pw.initPW(data_pol, th_mono, ph_mono);
				kinc[0] = pw.vW(0); kinc[1] = pw.vW(1); kinc[2] = pw.vW(2);
				einc[0] = pw.vE(0); einc[1] = pw.vE(1); einc[2] = pw.vE(2);
				hinc[0] = pw.vH(0); hinc[1] = pw.vH(1); hinc[2] = pw.vH(2);

				computeV(solver, kinc, einc, hinc);
				// 写入 V
				/*for (int i = 0; i < solver.row; ++i) {
					vFile << std::setw(16) << th_mono << "\t" << ph_mono << "\t" << i << "\t" << solver.Vm[i] << "\n";
				}
				vFile << "\n";*/

				// 求解
				int itr_max = solver.row, mr = 50;
				double tol_abs = 1.0e-8;
				double tol_rel = 1.0e-4;
				std::complex<double>* I_solve = new std::complex<double>[solver.row];
				std::complex<double>* Vec_R = new std::complex<double>[solver.row];
				for (int i = 0; i < solver.row; ++i) {
					I_solve[i] = 0.0;
					Vec_R[i] = solver.Vm[i];
				}
				solver.matrix_solver(solver.row, I_solve, Vec_R, itr_max, mr, tol_abs, tol_rel);
				// 写入 I
				/*for (int i = 0; i < solver.row; ++i) {
					iFile << std::setw(16) << th_mono << "\t" << ph_mono << "\t" << i << "\t" << I_solve[i] << "\n";
				}
				iFile << "\n";*/

				// 计算 RCS
				std::vector<std::complex<double>> RCS_pre(3, 0.0);
				RCSUtils::calculateRCS(solver, I_solve, RCS_pre, solver.wave.k1(), solver.wave.eta1(), th_mono, ph_mono);
				double vec_v[3] = { pw.vTh(0), pw.vTh(1), pw.vTh(2) };
				double vec_h[3] = { pw.vPh(0), pw.vPh(1), pw.vPh(2) };
				std::complex<double> E_v = vec_v[0] * RCS_pre[0] + vec_v[1] * RCS_pre[1] + vec_v[2] * RCS_pre[2];
				std::complex<double> E_h = vec_h[0] * RCS_pre[0] + vec_h[1] * RCS_pre[1] + vec_h[2] * RCS_pre[2];
				double rcs = 10.0 * log10((norm(E_v) + norm(E_h)) / (4.0 * Pi));
				double rcs_v = 10.0 * log10(norm(E_v) / (4.0 * Pi));
				double rcs_h = 10.0 * log10(norm(E_h) / (4.0 * Pi));
				rcsFile << std::fixed << std::setprecision(12);
				rcsFile << std::setw(16) << th_mono << "\t" << ph_mono << "\t" << idx << "\t" << rcs << "\t" << rcs_v << "\t" << rcs_h << std::endl;
				std::cout << "Mono angle point " << (idx + 1) << " / " << totalPoints
					<< " completed: theta=" << th_mono << ", phi=" << ph_mono << std::endl;

				delete[] I_solve;
				delete[] Vec_R;
			}
		}
		rcsFile.close();
		/*iFile.close();
		vFile.close();*/
		std::cout << "RCS export completed." << std::endl;
	}

	template<typename SolverType, typename ComputeVFunc>
	inline void computeDualStatic_PMCHWT(SolverType& solver, const RCSExportConfig& cfg, ComputeVFunc computeV) {
		double inc_th = cfg.incTheta;
		double inc_ph = cfg.incPhi;
		std::string pol_wave = cfg.polWaveStr();
		double data_pol = (pol_wave == "V") ? 0.0 : 90.0;

		double kinc[3], einc[3], hinc[3];
		EMSource pw;
		pw.initPW(data_pol, inc_th, inc_ph);
		kinc[0] = pw.vW(0); kinc[1] = pw.vW(1); kinc[2] = pw.vW(2);
		einc[0] = pw.vE(0); einc[1] = pw.vE(1); einc[2] = pw.vE(2);
		hinc[0] = pw.vH(0); hinc[1] = pw.vH(1); hinc[2] = pw.vH(2);

		computeV(solver, kinc, einc, hinc);

		//exportV(solver, cfg);

		// PMCHW: 矩阵大小 = 2N = 2 * solver.row
		int totalRow = 2 * solver.row;
		int itr_max = 2 * solver.row, mr = 50;
		double tol_abs = 1.0e-8;
		double tol_rel = 1.0e-4;
		std::complex<double>* I_solve = new std::complex<double>[totalRow];
		std::complex<double>* Vec_R = new std::complex<double>[totalRow];
		for (int i = 0; i < totalRow; ++i) {
			I_solve[i] = std::complex<double>(0.0, 0.0);
			Vec_R[i] = solver.Vm[i];
		}
		solver.matrix_solver(totalRow, I_solve, Vec_R,
			itr_max, mr, tol_abs, tol_rel);

		// 前 N = I_J, 后 N = I_M
		std::complex<double>* I_J = I_solve;
		std::complex<double>* I_M = I_solve + solver.row;

		auto rcsFile = cfg.openFile("_RCS", "txt");
		if (!rcsFile.is_open()) throw std::runtime_error("Cannot open RCS file");
		const int thetaNum = calculateAngleCount(cfg.scaThetaStart, cfg.scaThetaEnd, cfg.scaStep);
		const int phiNum = calculateAngleCount(cfg.scaPhiStart, cfg.scaPhiEnd, cfg.scaStep);
		const int totalPoints = thetaNum * phiNum;
		int idx = 0;
		for (int it = 0; it < thetaNum; ++it) {
			double th_dual = calculateAngleAt(
				cfg.scaThetaStart, cfg.scaThetaEnd, cfg.scaStep, it, thetaNum);
			for (int ip = 0; ip < phiNum; ++ip, ++idx) {
				double ph_dual = calculateAngleAt(
					cfg.scaPhiStart, cfg.scaPhiEnd, cfg.scaStep, ip, phiNum);
				std::vector<std::complex<double>> RCS_pre(3, 0.0);
				RCSUtils::calculateRCS_PMCHWT(solver, I_J, I_M, RCS_pre,
					solver.wave.k1(), solver.wave.k2(), solver.wave.eta1(), solver.wave.eta2(),
					th_dual, ph_dual);

				pw.initPW(data_pol, th_dual, ph_dual);
				double vec_v[3] = { pw.vTh(0), pw.vTh(1), pw.vTh(2) };
				double vec_h[3] = { pw.vPh(0), pw.vPh(1), pw.vPh(2) };
				std::complex<double> E_v = vec_v[0] * RCS_pre[0]
					+ vec_v[1] * RCS_pre[1] + vec_v[2] * RCS_pre[2];
				std::complex<double> E_h = vec_h[0] * RCS_pre[0]
					+ vec_h[1] * RCS_pre[1] + vec_h[2] * RCS_pre[2];
				double rcs = 10.0 * log10((norm(E_v) + norm(E_h)) / (4.0 * Pi));
				double rcs_v = 10.0 * log10(norm(E_v) / (4.0 * Pi));
				double rcs_h = 10.0 * log10(norm(E_h) / (4.0 * Pi));
				rcsFile << std::fixed << std::setprecision(12);
				rcsFile << std::setw(16) << th_dual << "\t" << ph_dual << "\t"
					<< idx << "\t" << rcs << "\t" << rcs_v << "\t" << rcs_h
					<< std::endl;
			}
		}
		size_t mem = solver.computeMem();
		delete[] I_solve;
		delete[] Vec_R;
		rcsFile << "\n\nMoM memory = "
			<< mem / (1024.0 * 1024.0) << " MB\n";
		rcsFile.close();
		std::cout << "MoM memory = "
			<< mem / (1024.0 * 1024.0) << " MB" << std::endl;
		std::cout << "RCS (PMCHW Dual) export completed." << std::endl;
	}

	template<typename SolverType, typename ComputeVFunc>
	inline void computeMonoStatic_PMCHWT(SolverType& solver, const RCSExportConfig& cfg, ComputeVFunc computeV) {
		auto rcsFile = cfg.openFile("_RCS", "txt");
		if (!rcsFile.is_open()) throw std::runtime_error("Cannot open RCS file");
		const int thetaNum = calculateAngleCount(cfg.scaThetaStart, cfg.scaThetaEnd, cfg.scaStep);
		const int phiNum = calculateAngleCount(cfg.scaPhiStart, cfg.scaPhiEnd, cfg.scaStep);
		const int totalPoints = thetaNum * phiNum;
		int totalRow = 2 * solver.row;
		int idx = 0;
		for (int it = 0; it < thetaNum; ++it) {
			double th_mono = calculateAngleAt(
				cfg.scaThetaStart, cfg.scaThetaEnd, cfg.scaStep, it, thetaNum);
			for (int ip = 0; ip < phiNum; ++ip, ++idx) {
				double ph_mono = calculateAngleAt(
					cfg.scaPhiStart, cfg.scaPhiEnd, cfg.scaStep, ip, phiNum);
				std::string pol_wave = cfg.polWaveStr();
				double data_pol = (pol_wave == "V") ? 0.0 : 90.0;

				double kinc[3], einc[3], hinc[3];
				EMSource pw;
				pw.initPW(data_pol, th_mono, ph_mono);
				kinc[0] = pw.vW(0); kinc[1] = pw.vW(1); kinc[2] = pw.vW(2);
				einc[0] = pw.vE(0); einc[1] = pw.vE(1); einc[2] = pw.vE(2);
				hinc[0] = pw.vH(0); hinc[1] = pw.vH(1); hinc[2] = pw.vH(2);

				computeV(solver, kinc, einc, hinc);

				int itr_max = 2 * solver.row, mr = 50;
				double tol_abs = 1.0e-8;
				double tol_rel = 1.0e-4;
				std::complex<double>* I_solve = new std::complex<double>[totalRow];
				std::complex<double>* Vec_R = new std::complex<double>[totalRow];
				for (int i = 0; i < totalRow; ++i) {
					I_solve[i] = 0.0;
					Vec_R[i] = solver.Vm[i];
				}
				solver.matrix_solver(totalRow, I_solve, Vec_R,
					itr_max, mr, tol_abs, tol_rel);

				std::complex<double>* I_J = I_solve;
				std::complex<double>* I_M = I_solve + solver.row;

				std::vector<std::complex<double>> RCS_pre(3, 0.0);
				RCSUtils::calculateRCS_PMCHWT(solver, I_J, I_M, RCS_pre,
					solver.wave.k1(), solver.wave.k2(), solver.wave.eta1(), solver.wave.eta2(),
					th_mono, ph_mono);

				double vec_v[3] = { pw.vTh(0), pw.vTh(1), pw.vTh(2) };
				double vec_h[3] = { pw.vPh(0), pw.vPh(1), pw.vPh(2) };
				std::complex<double> E_v = vec_v[0] * RCS_pre[0]
					+ vec_v[1] * RCS_pre[1] + vec_v[2] * RCS_pre[2];
				std::complex<double> E_h = vec_h[0] * RCS_pre[0]
					+ vec_h[1] * RCS_pre[1] + vec_h[2] * RCS_pre[2];
				double rcs = 10.0 * log10((norm(E_v) + norm(E_h)) / (4.0 * Pi));
				double rcs_v = 10.0 * log10(norm(E_v) / (4.0 * Pi));
				double rcs_h = 10.0 * log10(norm(E_h) / (4.0 * Pi));
				rcsFile << std::fixed << std::setprecision(12);
				rcsFile << std::setw(16) << th_mono << "\t" << ph_mono << "\t"
					<< idx << "\t" << rcs << "\t" << rcs_v << "\t" << rcs_h
					<< std::endl;
				std::cout << "PMCHWT mono angle point " << (idx + 1) << " / " << totalPoints
					<< " completed: theta=" << th_mono << ", phi=" << ph_mono << std::endl;

				delete[] I_solve;
				delete[] Vec_R;
			}
		}
		rcsFile.close();
		std::cout << "RCS (PMCHW Mono) export completed." << std::endl;
	}

	template<typename SolverType>
	inline void exportHSBSamples(SolverType& solver, const RCSExportConfig& cfg,
		const std::vector<std::vector<HSBMonoSample>>& rings)
	{
		(void)solver;
		auto sampleFile = cfg.openFile("_FIELD_HSB_SAMPLES", "txt");
		if (!sampleFile.is_open()) {
			throw std::runtime_error("Cannot open HSB samples field file");
		}

		const std::string pol_wave = cfg.polWaveStr();
		const double data_pol = (pol_wave == "V") ? 0.0 : 90.0;
		sampleFile << "sample_idx\tring_idx\tring_sample_idx\ttheta_deg\tphi_deg\tphi_order\t"
			<< "Ev_re\tEv_im\tEh_re\tEh_im\trcs_db\trcs_v_db\trcs_h_db\n";
		sampleFile << std::fixed << std::setprecision(12);

		int sampleIdx = 0;
		for (int ringIdx = 0; ringIdx < static_cast<int>(rings.size()); ++ringIdx) {
			const auto& ring = rings[ringIdx];
			for (int ringSampleIdx = 0;
				ringSampleIdx < static_cast<int>(ring.size());
				++ringSampleIdx) {
				const HSBMonoSample& sample = ring[ringSampleIdx];
				const std::vector<std::complex<double>>& RCS_pre = sample.field;

				EMSource pw;
				pw.initPW(data_pol, sample.thetaDeg, sample.phiDeg);
				double vec_v[3] = { pw.vTh(0), pw.vTh(1), pw.vTh(2) };
				double vec_h[3] = { pw.vPh(0), pw.vPh(1), pw.vPh(2) };

				std::complex<double> E_v =
					vec_v[0] * RCS_pre[0] +
					vec_v[1] * RCS_pre[1] +
					vec_v[2] * RCS_pre[2];
				std::complex<double> E_h =
					vec_h[0] * RCS_pre[0] +
					vec_h[1] * RCS_pre[1] +
					vec_h[2] * RCS_pre[2];

				double rcs = 10.0 * std::log10((std::norm(E_v) + std::norm(E_h)) / (4.0 * Pi));
				double rcs_v = 10.0 * std::log10(std::norm(E_v) / (4.0 * Pi));
				double rcs_h = 10.0 * std::log10(std::norm(E_h) / (4.0 * Pi));

				sampleFile << sampleIdx << "\t"
					<< ringIdx << "\t" << ringSampleIdx << "\t"
					<< sample.thetaDeg << "\t" << sample.phiDeg << "\t"
					<< sample.phiOrder << "\t"
					<< E_v.real() << "\t" << E_v.imag() << "\t"
					<< E_h.real() << "\t" << E_h.imag() << "\t"
					<< rcs << "\t" << rcs_v << "\t" << rcs_h << "\n";
				++sampleIdx;
			}
		}
	}

	template<typename SolverType>
	inline void exportHSBScanGrid(SolverType& solver, const RCSExportConfig& cfg,
		const std::vector<std::vector<HSBMonoSample>>& rings)
	{
		auto angleCount = [](double start, double end, double step) {
			const double diff = std::abs(end - start);
			if (diff < 1.0e-9) {
				return 1;
			}
			return static_cast<int>(std::ceil(diff / step - 1.0e-9)) + 1;
			};
		auto angleAt = [](double start, double end, double step, int idx, int count) {
			if (idx == count - 1) {
				return end;
			}
			const double direction = (end >= start) ? 1.0 : -1.0;
			return start + direction * idx * step;
			};

		(void)solver;
		auto rcsFile = cfg.openFile("_RCS_HSB", "txt");
		if (!rcsFile.is_open()) {
			throw std::runtime_error("Cannot open HSB RCS file");
		}

		std::ofstream fieldFile;
		if (cfg.exportHSBComplexField) {
			fieldFile = cfg.openFile("_FIELD_HSB", "txt");
			if (!fieldFile.is_open()) {
				throw std::runtime_error("Cannot open HSB complex field file");
			}
		}

		const std::string pol_wave = cfg.polWaveStr();
		const double data_pol = (pol_wave == "V") ? 0.0 : 90.0;
		const int thetaNum = angleCount(cfg.scaThetaStart, cfg.scaThetaEnd, cfg.scaStep);
		const int phiNum = angleCount(cfg.scaPhiStart, cfg.scaPhiEnd, cfg.scaStep);

		rcsFile << std::fixed << std::setprecision(12);
		if (cfg.exportHSBComplexField) {
			fieldFile << std::fixed << std::setprecision(12);
		}

		for (int it = 0; it < thetaNum; ++it) {
			double thetaDeg = angleAt(
				cfg.scaThetaStart, cfg.scaThetaEnd, cfg.scaStep, it, thetaNum);

			for (int ip = 0; ip < phiNum; ++ip) {
				double phiDeg = angleAt(
					cfg.scaPhiStart, cfg.scaPhiEnd, cfg.scaStep, ip, phiNum);
				const double phiInterpDeg = my_hsb::normalizePhiDeg(phiDeg);

				std::vector<std::complex<double>> RCS_pre =
					my_hsb::interpolateHSBField(rings, thetaDeg, phiInterpDeg);

				EMSource pw;
				pw.initPW(data_pol, thetaDeg, phiDeg);
				double vec_v[3] = { pw.vTh(0), pw.vTh(1), pw.vTh(2) };
				double vec_h[3] = { pw.vPh(0), pw.vPh(1), pw.vPh(2) };

				std::complex<double> E_v =
					vec_v[0] * RCS_pre[0] +
					vec_v[1] * RCS_pre[1] +
					vec_v[2] * RCS_pre[2];
				std::complex<double> E_h =
					vec_h[0] * RCS_pre[0] +
					vec_h[1] * RCS_pre[1] +
					vec_h[2] * RCS_pre[2];

				double rcs = 10.0 * std::log10((std::norm(E_v) + std::norm(E_h)) / (4.0 * Pi));
				double rcs_v = 10.0 * std::log10(std::norm(E_v) / (4.0 * Pi));
				double rcs_h = 10.0 * std::log10(std::norm(E_h) / (4.0 * Pi));

				rcsFile << std::setw(16) << thetaDeg << "\t" << phiDeg << "\t"
					<< it << "\t" << ip << "\t"
					<< rcs << "\t" << rcs_v << "\t" << rcs_h << "\n";

				if (cfg.exportHSBComplexField) {
					fieldFile << std::setw(16) << thetaDeg << "\t" << phiDeg << "\t"
						<< it << "\t" << ip << "\t"
						<< E_v.real() << "\t" << E_v.imag() << "\t"
						<< E_h.real() << "\t" << E_h.imag() << "\n";
				}
			}
		}
	}
	template<typename SolverType, typename ComputeVFunc>
	inline void computeMonoStatic_HSB(SolverType& solver, const RCSExportConfig& cfg,
		ComputeVFunc computeV, double rho = 1.2)
	{
		const std::string pol_wave = cfg.polWaveStr();
		const double data_pol = (pol_wave == "V") ? 0.0 : 90.0;
		const bool isPMCHWT = (solver.integralEquType_ == 2);
		const int totalRow = isPMCHWT ? 2 * solver.row : solver.row;
		const double radius = my_hsb::estimateBoundingSphereRadius(solver);
		const double W = std::max(1.0, rho * solver.wave.k1_abs() * radius);
		const int thetaCount = std::max(2, static_cast<int>(std::ceil(2.0 * W)));
		int expectedSampleCount = 0;
		for (int m = 1; m <= thetaCount; ++m) {
			const double theta = Pi * static_cast<double>(m) /
				static_cast<double>(thetaCount + 1);
			const int phiOrder = std::max(
				1, static_cast<int>(std::ceil(W * std::sin(theta))));
			expectedSampleCount += 2 * phiOrder + 1;
		}
		std::cout << "HSB accurate sample points to compute = "
			<< expectedSampleCount << std::endl;

		std::vector<std::vector<HSBMonoSample>> rings;
		rings.reserve(thetaCount);

		int sampleCount = 0;
		for (int m = 1; m <= thetaCount; ++m) {
			const double theta = Pi * static_cast<double>(m) /
				static_cast<double>(thetaCount + 1);
			const double thetaDeg = theta * 180.0 / Pi;
			const int phiOrder = std::max(
				1, static_cast<int>(std::ceil(W * std::sin(theta))));
			const int phiCount = 2 * phiOrder + 1;

			std::vector<HSBMonoSample> ring;
			ring.reserve(phiCount);

			for (int n = -phiOrder; n <= phiOrder; ++n) {
				double phiDeg = my_hsb::normalizePhiDeg(
					360.0 * static_cast<double>(n) /
					static_cast<double>(phiCount));

				EMSource pw;
				pw.initPW(data_pol, thetaDeg, phiDeg);
				double kinc[3] = { pw.vW(0), pw.vW(1), pw.vW(2) };
				double einc[3] = { pw.vE(0), pw.vE(1), pw.vE(2) };
				double hinc[3] = { pw.vH(0), pw.vH(1), pw.vH(2) };

				computeV(solver, kinc, einc, hinc);

				std::vector<std::complex<double>> I_solve(
					totalRow, std::complex<double>(0.0, 0.0));
				std::vector<std::complex<double>> Vec_R(totalRow);
				for (int i = 0; i < totalRow; ++i) {
					Vec_R[i] = solver.Vm[i];
				}

				int itr_max = totalRow;
				int mr = 50;
				double tol_abs = 1.0e-8;
				double tol_rel = 1.0e-4;
				solver.matrix_solver(
					totalRow, I_solve.data(), Vec_R.data(),
					itr_max, mr, tol_abs, tol_rel);
				std::cout << "HSB accurate sample point "
					<< (sampleCount + 1) << " / " << expectedSampleCount
					<< " solved" << std::endl;

				std::vector<std::complex<double>> RCS_pre(3, 0.0);
				if (isPMCHWT) {
					RCSUtils::calculateRCS_PMCHWT(
						solver, I_solve.data(), I_solve.data() + solver.row,
						RCS_pre, solver.wave.k1(), solver.wave.k2(),
						solver.wave.eta1(), solver.wave.eta2(),
						thetaDeg, phiDeg);
				}
				else {
					RCSUtils::calculateRCS(
						solver, I_solve.data(), RCS_pre,
						solver.wave.k1(), solver.wave.eta1(),
						thetaDeg, phiDeg);
				}

				HSBMonoSample sample;
				sample.thetaDeg = thetaDeg;
				sample.phiDeg = phiDeg;
				sample.phiOrder = phiOrder;
				sample.field = std::move(RCS_pre);
				ring.push_back(std::move(sample));
				++sampleCount;
			}

			rings.push_back(std::move(ring));
		}

		if (cfg.exportHSBSamples) {
			RCSUtils::exportHSBSamples(solver, cfg, rings);
		}
		RCSUtils::exportHSBScanGrid(solver, cfg, rings);

		const size_t mem = solver.computeMem();
		std::cout << "HSB sample points = " << sampleCount << std::endl;
		std::cout << "Algorithm memory = "
			<< mem / (1024.0 * 1024.0) << " MB" << std::endl;
		std::cout << "RCS HSB export completed." << std::endl;
	}
}
