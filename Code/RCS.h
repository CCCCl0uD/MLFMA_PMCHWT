#pragma once

#include "RCSExportConfig.h"
#include "EMSource.h"
#include "MovingTarget.h"
#include "ACARHS.h"

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <utility>
#include <vector>

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
	struct RCSGridValue {
		double rcs = 0.0;
		double rcs_v = 0.0;
		double rcs_h = 0.0;
		std::complex<double> E_v = 0.0;
		std::complex<double> E_h = 0.0;
	};

	inline bool isThetaPole(double thetaDeg) {
		const double eps = 1.0e-9;
		return std::abs(thetaDeg) < eps || std::abs(thetaDeg - 180.0) < eps;
	}

	inline double canonicalPolePhiDeg() {
		return 0.0;
	}

	inline RCSGridValue projectRCSGridValue(
		const std::vector<std::complex<double>>& RCS_pre,
		double data_pol, double thetaDeg, double phiDeg)
	{
		EMSource pw;
		pw.initPW(data_pol, thetaDeg, phiDeg);
		double vec_v[3] = { pw.vTh(0), pw.vTh(1), pw.vTh(2) };
		double vec_h[3] = { pw.vPh(0), pw.vPh(1), pw.vPh(2) };

		RCSGridValue value;
		value.E_v = vec_v[0] * RCS_pre[0] + vec_v[1] * RCS_pre[1] + vec_v[2] * RCS_pre[2];
		value.E_h = vec_h[0] * RCS_pre[0] + vec_h[1] * RCS_pre[1] + vec_h[2] * RCS_pre[2];
		value.rcs = 10.0 * std::log10((std::norm(value.E_v) + std::norm(value.E_h)) / (4.0 * Pi));
		value.rcs_v = 10.0 * std::log10(std::norm(value.E_v) / (4.0 * Pi));
		value.rcs_h = 10.0 * std::log10(std::norm(value.E_h) / (4.0 * Pi));
		return value;
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

				// ���湱��
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

				// ���湱��
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
		// �� cfg ��ȡ����
		double inc_th = cfg.incTheta;
		double inc_ph = cfg.incPhi;
		std::string pol_wave = cfg.polWaveStr(); // "V" or "H"
		double data_pol = (pol_wave == "V") ? 0.0 : 90.0;

		// �������䲨
		double kinc[3], einc[3], hinc[3];
		EMSource pw;
		pw.initPW(data_pol, inc_th, inc_ph);
		kinc[0] = pw.vW(0); kinc[1] = pw.vW(1); kinc[2] = pw.vW(2);
		einc[0] = pw.vE(0); einc[1] = pw.vE(1); einc[2] = pw.vE(2);
		hinc[0] = pw.vH(0); hinc[1] = pw.vH(1); hinc[2] = pw.vH(2);

		// �����ұ�����
		computeV(solver, kinc, einc, hinc);

		// ���� V
		//exportV(solver, cfg);

		// ���
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

		// ���� I
		//exportI(solver, cfg, I_solve);

		// ����˫վ RCS
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

				// ���㼫������
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
			const bool thetaPole = isThetaPole(th_mono);
			bool hasPoleValue = false;
			RCSGridValue poleValue;
			for (int ip = 0; ip < phiNum; ++ip, ++idx) {
				double ph_mono = calculateAngleAt(
					cfg.scaPhiStart, cfg.scaPhiEnd, cfg.scaStep, ip, phiNum);
				if (thetaPole && hasPoleValue) {
					rcsFile << std::fixed << std::setprecision(12);
					rcsFile << std::setw(16) << th_mono << "\t" << ph_mono << "\t" << idx << "\t"
						<< poleValue.rcs << "\t" << poleValue.rcs_v << "\t" << poleValue.rcs_h << std::endl;
					std::cout << "Mono pole angle point " << (idx + 1) << " / " << totalPoints
						<< " copied: theta=" << th_mono << ", phi=" << ph_mono << std::endl;
					continue;
				}
				const double ph_compute = thetaPole ? canonicalPolePhiDeg() : ph_mono;
				std::string pol_wave = cfg.polWaveStr();
				double data_pol = (pol_wave == "V") ? 0.0 : 90.0;

				double kinc[3], einc[3], hinc[3];
				EMSource pw;
				pw.initPW(data_pol, th_mono, ph_compute);
				kinc[0] = pw.vW(0); kinc[1] = pw.vW(1); kinc[2] = pw.vW(2);
				einc[0] = pw.vE(0); einc[1] = pw.vE(1); einc[2] = pw.vE(2);
				hinc[0] = pw.vH(0); hinc[1] = pw.vH(1); hinc[2] = pw.vH(2);

				computeV(solver, kinc, einc, hinc);
				// д�� V
				/*for (int i = 0; i < solver.row; ++i) {
					vFile << std::setw(16) << th_mono << "\t" << ph_mono << "\t" << i << "\t" << solver.Vm[i] << "\n";
				}
				vFile << "\n";*/

				// ���
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
				// д�� I
				/*for (int i = 0; i < solver.row; ++i) {
					iFile << std::setw(16) << th_mono << "\t" << ph_mono << "\t" << i << "\t" << I_solve[i] << "\n";
				}
				iFile << "\n";*/

				// ���� RCS
				std::vector<std::complex<double>> RCS_pre(3, 0.0);
				RCSUtils::calculateRCS(solver, I_solve, RCS_pre, solver.wave.k1(), solver.wave.eta1(), th_mono, ph_compute);
				double vec_v[3] = { pw.vTh(0), pw.vTh(1), pw.vTh(2) };
				double vec_h[3] = { pw.vPh(0), pw.vPh(1), pw.vPh(2) };
				std::complex<double> E_v = vec_v[0] * RCS_pre[0] + vec_v[1] * RCS_pre[1] + vec_v[2] * RCS_pre[2];
				std::complex<double> E_h = vec_h[0] * RCS_pre[0] + vec_h[1] * RCS_pre[1] + vec_h[2] * RCS_pre[2];
				double rcs = 10.0 * log10((norm(E_v) + norm(E_h)) / (4.0 * Pi));
				double rcs_v = 10.0 * log10(norm(E_v) / (4.0 * Pi));
				double rcs_h = 10.0 * log10(norm(E_h) / (4.0 * Pi));
				if (thetaPole) {
					poleValue.rcs = rcs;
					poleValue.rcs_v = rcs_v;
					poleValue.rcs_h = rcs_h;
					hasPoleValue = true;
				}
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

		// PMCHW: �����С = 2N = 2 * solver.row
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

		// ǰ N = I_J, �� N = I_M
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
			const bool thetaPole = isThetaPole(th_mono);
			bool hasPoleValue = false;
			RCSGridValue poleValue;
			for (int ip = 0; ip < phiNum; ++ip, ++idx) {
				double ph_mono = calculateAngleAt(
					cfg.scaPhiStart, cfg.scaPhiEnd, cfg.scaStep, ip, phiNum);
				if (thetaPole && hasPoleValue) {
					rcsFile << std::fixed << std::setprecision(12);
					rcsFile << std::setw(16) << th_mono << "\t" << ph_mono << "\t" << idx << "\t"
						<< poleValue.rcs << "\t" << poleValue.rcs_v << "\t" << poleValue.rcs_h << std::endl;
					std::cout << "Mono pole angle point " << (idx + 1) << " / " << totalPoints
						<< " copied: theta=" << th_mono << ", phi=" << ph_mono << std::endl;
					continue;
				}
				const double ph_compute = thetaPole ? canonicalPolePhiDeg() : ph_mono;
				std::string pol_wave = cfg.polWaveStr();
				double data_pol = (pol_wave == "V") ? 0.0 : 90.0;

				double kinc[3], einc[3], hinc[3];
				EMSource pw;
				pw.initPW(data_pol, th_mono, ph_compute);
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
					th_mono, ph_compute);

				double vec_v[3] = { pw.vTh(0), pw.vTh(1), pw.vTh(2) };
				double vec_h[3] = { pw.vPh(0), pw.vPh(1), pw.vPh(2) };
				std::complex<double> E_v = vec_v[0] * RCS_pre[0]
					+ vec_v[1] * RCS_pre[1] + vec_v[2] * RCS_pre[2];
				std::complex<double> E_h = vec_h[0] * RCS_pre[0]
					+ vec_h[1] * RCS_pre[1] + vec_h[2] * RCS_pre[2];
				double rcs = 10.0 * log10((norm(E_v) + norm(E_h)) / (4.0 * Pi));
				double rcs_v = 10.0 * log10(norm(E_v) / (4.0 * Pi));
				double rcs_h = 10.0 * log10(norm(E_h) / (4.0 * Pi));
				if (thetaPole) {
					poleValue.rcs = rcs;
					poleValue.rcs_v = rcs_v;
					poleValue.rcs_h = rcs_h;
					hasPoleValue = true;
				}
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
	inline void exportFSBSamples(SolverType& solver, const RCSExportConfig& cfg,
		const std::vector<std::vector<FSBMonoSample>>& rings)
	{
		(void)solver;
		auto sampleFile = cfg.openFile("_FIELD_FSB_SAMPLES", "txt");
		if (!sampleFile.is_open()) {
			throw std::runtime_error("Cannot open FSB samples field file");
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
				const FSBMonoSample& sample = ring[ringSampleIdx];
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
	inline void exportFSBScanGrid(SolverType& solver, const RCSExportConfig& cfg,
		const std::vector<std::vector<FSBMonoSample>>& rings)
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
		auto rcsFile = cfg.openFile("_RCS_FSB", "txt");
		if (!rcsFile.is_open()) {
			throw std::runtime_error("Cannot open FSB RCS file");
		}

		std::ofstream fieldFile;
		if (cfg.exportFSBComplexField) {
			fieldFile = cfg.openFile("_FIELD_FSB", "txt");
			if (!fieldFile.is_open()) {
				throw std::runtime_error("Cannot open FSB complex field file");
			}
		}

		const std::string pol_wave = cfg.polWaveStr();
		const double data_pol = (pol_wave == "V") ? 0.0 : 90.0;
		const int thetaNum = angleCount(cfg.scaThetaStart, cfg.scaThetaEnd, cfg.scaStep);
		const int phiNum = angleCount(cfg.scaPhiStart, cfg.scaPhiEnd, cfg.scaStep);

		rcsFile << std::fixed << std::setprecision(12);
		if (cfg.exportFSBComplexField) {
			fieldFile << std::fixed << std::setprecision(12);
		}

		for (int it = 0; it < thetaNum; ++it) {
			double thetaDeg = angleAt(
				cfg.scaThetaStart, cfg.scaThetaEnd, cfg.scaStep, it, thetaNum);

			const bool thetaPole = isThetaPole(thetaDeg);
			bool hasPoleValue = false;
			RCSGridValue poleValue;

			for (int ip = 0; ip < phiNum; ++ip) {
				double phiDeg = angleAt(
					cfg.scaPhiStart, cfg.scaPhiEnd, cfg.scaStep, ip, phiNum);
				RCSGridValue value;
				if (thetaPole && hasPoleValue) {
					value = poleValue;
				}
				else {
					const double phiComputeDeg = thetaPole ? canonicalPolePhiDeg() : phiDeg;
					const double phiInterpDeg = my_fsb::normalizePhiDeg(phiComputeDeg);
					std::vector<std::complex<double>> RCS_pre =
						my_fsb::interpolateFSBField(rings, thetaDeg, phiInterpDeg);
					value = projectRCSGridValue(RCS_pre, data_pol, thetaDeg, phiComputeDeg);
					if (thetaPole) {
						poleValue = value;
						hasPoleValue = true;
					}
				}

				rcsFile << std::setw(16) << thetaDeg << "\t" << phiDeg << "\t"
					<< it << "\t" << ip << "\t"
					<< value.rcs << "\t" << value.rcs_v << "\t" << value.rcs_h << "\n";

				if (cfg.exportFSBComplexField) {
					fieldFile << std::setw(16) << thetaDeg << "\t" << phiDeg << "\t"
						<< it << "\t" << ip << "\t"
						<< value.E_v.real() << "\t" << value.E_v.imag() << "\t"
						<< value.E_h.real() << "\t" << value.E_h.imag() << "\n";
				}
			}
		}
	}
	template<typename SolverType>
	inline void exportMovingAngles(SolverType& solver, const RCSExportConfig& cfg,
		const std::vector<MovingAngleSample>& movingAngles)
	{
		(void)solver;
		auto angleFile = cfg.openFile("_ANGLES_FSB_MLFMA_MOVING", "txt");
		if (!angleFile.is_open()) {
			throw std::runtime_error("Cannot open FSB-MLFMA-Moving angle file");
		}
		angleFile << std::fixed << std::setprecision(12);
		angleFile << "idx\ttime_s\tpos_x\tpos_y\tpos_z\trange_m\t"
			<< "theta_global_deg\tphi_global_deg\ttheta_body_deg\tphi_body_deg\n";
		for (int i = 0; i < static_cast<int>(movingAngles.size()); ++i) {
			const auto& sample = movingAngles[i];
			angleFile << i << "\t" << sample.time_s << "\t"
				<< sample.pos[0] << "\t" << sample.pos[1] << "\t" << sample.pos[2] << "\t"
				<< sample.range << "\t"
				<< sample.thetaGlobalDeg << "\t" << sample.phiGlobalDeg << "\t"
				<< sample.thetaBodyDeg << "\t" << sample.phiBodyDeg << "\n";
		}
	}

	template<typename SolverType>
	inline void exportACARHSMovingAngles(SolverType& solver, const RCSExportConfig& cfg,
		const std::vector<MovingAngleSample>& movingAngles)
	{
		(void)solver;
		auto angleFile = cfg.openFile("_ANGLES_ACARHS_MLFMA_MOVING", "txt");
		if (!angleFile.is_open()) {
			throw std::runtime_error("Cannot open ACA-RHS-MLFMA-Moving angle file");
		}
		angleFile << std::fixed << std::setprecision(12);
		angleFile << "idx\ttime_s\tpos_x\tpos_y\tpos_z\trange_m\t"
			<< "theta_global_deg\tphi_global_deg\ttheta_body_deg\tphi_body_deg\n";
		for (int i = 0; i < static_cast<int>(movingAngles.size()); ++i) {
			const auto& sample = movingAngles[i];
			angleFile << i << "\t" << sample.time_s << "\t"
				<< sample.pos[0] << "\t" << sample.pos[1] << "\t" << sample.pos[2] << "\t"
				<< sample.range << "\t"
				<< sample.thetaGlobalDeg << "\t" << sample.phiGlobalDeg << "\t"
				<< sample.thetaBodyDeg << "\t" << sample.phiBodyDeg << "\n";
		}
	}

	template<typename SolverType>
	inline void exportAccurateMovingAngles(SolverType& solver, const RCSExportConfig& cfg,
		const std::vector<MovingAngleSample>& movingAngles)
	{
		(void)solver;
		auto angleFile = cfg.openFile("_ANGLES_ACCURATE_MLFMA_MOVING", "txt");
		if (!angleFile.is_open()) {
			throw std::runtime_error("Cannot open Accurate-MLFMA-Moving angle file");
		}
		angleFile << std::fixed << std::setprecision(12);
		angleFile << "idx\ttime_s\tpos_x\tpos_y\tpos_z\trange_m\t"
			<< "theta_global_deg\tphi_global_deg\ttheta_body_deg\tphi_body_deg\n";
		for (int i = 0; i < static_cast<int>(movingAngles.size()); ++i) {
			const auto& sample = movingAngles[i];
			angleFile << i << "\t" << sample.time_s << "\t"
				<< sample.pos[0] << "\t" << sample.pos[1] << "\t" << sample.pos[2] << "\t"
				<< sample.range << "\t"
				<< sample.thetaGlobalDeg << "\t" << sample.phiGlobalDeg << "\t"
				<< sample.thetaBodyDeg << "\t" << sample.phiBodyDeg << "\n";
		}
	}

	inline double complexVectorNorm2(const std::vector<std::complex<double>>& values)
	{
		double result = 0.0;
		for (const auto& value : values) {
			result += std::norm(value);
		}
		return result;
	}

	template<typename ComputeColumnFunc>
	inline void exportACARHSCheck(const RCSExportConfig& cfg,
		const std::vector<MovingAngleSample>& movingAngles,
		const ACARHS::Factors& factors,
		ComputeColumnFunc computeColumn)
	{
		if (!cfg.exportACARHSCheck || cfg.acaRHSCheckCount <= 0) {
			return;
		}

		auto checkFile = cfg.openFile("_ACA_RHS_CHECK", "txt");
		if (!checkFile.is_open()) {
			throw std::runtime_error("Cannot open ACA-RHS check file");
		}

		const int timeCount = static_cast<int>(movingAngles.size());
		const int checkCount = std::min(cfg.acaRHSCheckCount, timeCount);
		checkFile << std::fixed << std::setprecision(12);
		checkFile << "time_idx\ttime_s\trhs_exact_norm\trhs_error_norm\trelative_error\n";

		int lastTimeIndex = -1;
		for (int c = 0; c < checkCount; ++c) {
			int timeIndex = 0;
			if (checkCount == 1) {
				timeIndex = 0;
			}
			else {
				timeIndex = static_cast<int>(std::round(
					static_cast<double>(c) * static_cast<double>(timeCount - 1) /
					static_cast<double>(checkCount - 1)));
			}
			if (timeIndex == lastTimeIndex) {
				continue;
			}
			lastTimeIndex = timeIndex;

			std::vector<std::complex<double>> exact;
			computeColumn(timeIndex, exact);
			std::vector<std::complex<double>> approx =
				ACARHS::reconstructColumn(factors, timeIndex);
			if (exact.size() != approx.size()) {
				throw std::runtime_error("ACA-RHS check vector size mismatch.");
			}

			std::vector<std::complex<double>> error(exact.size());
			for (int i = 0; i < static_cast<int>(exact.size()); ++i) {
				error[i] = exact[i] - approx[i];
			}

			const double exactNorm = std::sqrt(complexVectorNorm2(exact));
			const double errorNorm = std::sqrt(complexVectorNorm2(error));
			const double relativeError = (exactNorm > 0.0) ? errorNorm / exactNorm : errorNorm;

			checkFile << timeIndex << "\t"
				<< movingAngles[timeIndex].time_s << "\t"
				<< exactNorm << "\t"
				<< errorNorm << "\t"
				<< relativeError << "\n";
		}
	}

	inline void exportACARHSFactors(const RCSExportConfig& cfg,
		const ACARHS::Factors& factors)
	{
		if (!cfg.exportACARHSFactors) {
			return;
		}

		auto uFile = cfg.openFile("_ACA_RHS_FACTORS_U", "txt");
		if (!uFile.is_open()) {
			throw std::runtime_error("Cannot open ACA-RHS U factor file");
		}
		uFile << std::fixed << std::setprecision(17);
		uFile << "rank_idx\trow_idx\treal\timag\n";
		for (int m = 0; m < factors.rank; ++m) {
			for (int row = 0; row < factors.rows; ++row) {
				const auto value = factors.U[static_cast<size_t>(m) * factors.rows + row];
				uFile << m << "\t" << row << "\t"
					<< value.real() << "\t" << value.imag() << "\n";
			}
		}

		auto vFile = cfg.openFile("_ACA_RHS_FACTORS_V", "txt");
		if (!vFile.is_open()) {
			throw std::runtime_error("Cannot open ACA-RHS V factor file");
		}
		vFile << std::fixed << std::setprecision(17);
		vFile << "rank_idx\ttime_idx\treal\timag\n";
		for (int m = 0; m < factors.rank; ++m) {
			for (int col = 0; col < factors.cols; ++col) {
				const auto value = factors.V[static_cast<size_t>(m) * factors.cols + col];
				vFile << m << "\t" << col << "\t"
					<< value.real() << "\t" << value.imag() << "\n";
			}
		}
	}

	template<typename SolverType, typename ComputeVFunc>
	inline void computeMonoStatic_Accurate_Moving(SolverType& solver, const RCSExportConfig& cfg,
		ComputeVFunc computeV)
	{
		if (solver.integralEquType_ == 2) {
			throw std::runtime_error("Accurate-MLFMA-Moving supports PEC mono only in this version.");
		}

		const std::string pol_wave = cfg.polWaveStr();
		const double data_pol = (pol_wave == "V") ? 0.0 : 90.0;
		const std::array<double, 3> radarPos = { cfg.radarX, cfg.radarY, cfg.radarZ };

		std::vector<MovingState> states = MovingTarget::readMovingPathFile(cfg.movingPathFile);
		std::vector<MovingAngleSample> movingAngles =
			MovingTarget::computeMovingAngleSamples(states, radarPos);
		if (movingAngles.empty()) {
			throw std::runtime_error("Accurate-MLFMA-Moving has no moving states.");
		}
		if (cfg.exportMovingAngles) {
			RCSUtils::exportAccurateMovingAngles(solver, cfg, movingAngles);
		}

		const int timeCount = static_cast<int>(movingAngles.size());
		const int rowCount = solver.row;
		std::cout << "Accurate-MLFMA-Moving states = " << timeCount << std::endl;
		std::cout << "Accurate-MLFMA-Moving RHS rows = " << rowCount << std::endl;

		auto rcsFile = cfg.openFile("_RCS_ACCURATE_MLFMA_MOVING", "txt");
		if (!rcsFile.is_open()) {
			throw std::runtime_error("Cannot open Accurate-MLFMA-Moving RCS file");
		}
		rcsFile << std::fixed << std::setprecision(12);
		rcsFile << "idx\ttime_s\trange_m\ttheta_global_deg\tphi_global_deg\t"
			<< "theta_body_deg\tphi_body_deg\trcs_db\trcs_v_db\trcs_h_db\n";

		std::ofstream fieldFile;
		if (cfg.exportMovingComplexField) {
			fieldFile = cfg.openFile("_FIELD_ACCURATE_MLFMA_MOVING", "txt");
			if (!fieldFile.is_open()) {
				throw std::runtime_error("Cannot open Accurate-MLFMA-Moving field file");
			}
			fieldFile << std::fixed << std::setprecision(12);
			fieldFile << "idx\ttime_s\trange_m\ttheta_body_deg\tphi_body_deg\t"
				<< "Ev_re\tEv_im\tEh_re\tEh_im\n";
		}

		std::vector<std::complex<double>> I_solve(
			rowCount, std::complex<double>(0.0, 0.0));
		std::vector<std::complex<double>> Vec_R(rowCount);

		for (int timeIndex = 0; timeIndex < timeCount; ++timeIndex) {
			const auto& angle = movingAngles[timeIndex];
			const double phiComputeDeg = isThetaPole(angle.thetaBodyDeg) ?
				canonicalPolePhiDeg() : angle.phiBodyDeg;

			EMSource pw;
			pw.initPW(data_pol, angle.thetaBodyDeg, phiComputeDeg);
			double kinc[3] = { pw.vW(0), pw.vW(1), pw.vW(2) };
			double einc[3] = { pw.vE(0), pw.vE(1), pw.vE(2) };
			double hinc[3] = { pw.vH(0), pw.vH(1), pw.vH(2) };

			computeV(solver, kinc, einc, hinc);
			std::fill(I_solve.begin(), I_solve.end(), std::complex<double>(0.0, 0.0));
			for (int row = 0; row < rowCount; ++row) {
				Vec_R[row] = solver.Vm[row];
			}

			int itr_max = rowCount;
			int mr = 50;
			double tol_abs = 1.0e-8;
			double tol_rel = 1.0e-4;
			solver.matrix_solver(
				rowCount, I_solve.data(), Vec_R.data(),
				itr_max, mr, tol_abs, tol_rel);

			std::vector<std::complex<double>> RCS_pre(3, 0.0);
			RCSUtils::calculateRCS(
				solver, I_solve.data(), RCS_pre,
				solver.wave.k1(), solver.wave.eta1(),
				angle.thetaBodyDeg, phiComputeDeg);
			RCSGridValue value = projectRCSGridValue(
				RCS_pre, data_pol, angle.thetaBodyDeg, phiComputeDeg);

			rcsFile << timeIndex << "\t" << angle.time_s << "\t" << angle.range << "\t"
				<< angle.thetaGlobalDeg << "\t" << angle.phiGlobalDeg << "\t"
				<< angle.thetaBodyDeg << "\t" << angle.phiBodyDeg << "\t"
				<< value.rcs << "\t" << value.rcs_v << "\t" << value.rcs_h << "\n";

			if (cfg.exportMovingComplexField) {
				fieldFile << timeIndex << "\t" << angle.time_s << "\t" << angle.range << "\t"
					<< angle.thetaBodyDeg << "\t" << angle.phiBodyDeg << "\t"
					<< value.E_v.real() << "\t" << value.E_v.imag() << "\t"
					<< value.E_h.real() << "\t" << value.E_h.imag() << "\n";
			}

			std::cout << "Accurate-MLFMA-Moving time point "
				<< (timeIndex + 1) << " / " << timeCount << " RCS computed" << std::endl;
		}

		const size_t mem = solver.computeMem() +
			my_memory::memory1D<std::complex<double>>(I_solve) +
			my_memory::memory1D<std::complex<double>>(Vec_R);
		std::cout << "Algorithm memory = "
			<< mem / (1024.0 * 1024.0) << " MB" << std::endl;
		std::cout << "RCS Accurate-MLFMA-Moving export completed." << std::endl;
	}

	template<typename SolverType, typename ComputeVFunc, typename ComputeVEntryFunc>
	inline void computeMonoStatic_ACARHS_Moving(SolverType& solver, const RCSExportConfig& cfg,
		ComputeVFunc computeV, ComputeVEntryFunc computeVEntry)
	{
		if (solver.integralEquType_ == 2) {
			throw std::runtime_error("ACA-RHS-MLFMA-Moving supports PEC mono only in this version.");
		}

		const std::string pol_wave = cfg.polWaveStr();
		const double data_pol = (pol_wave == "V") ? 0.0 : 90.0;
		const std::array<double, 3> radarPos = { cfg.radarX, cfg.radarY, cfg.radarZ };

		std::vector<MovingState> states = MovingTarget::readMovingPathFile(cfg.movingPathFile);
		std::vector<MovingAngleSample> movingAngles =
			MovingTarget::computeMovingAngleSamples(states, radarPos);
		if (movingAngles.empty()) {
			throw std::runtime_error("ACA-RHS-MLFMA-Moving has no moving states.");
		}
		if (cfg.exportMovingAngles) {
			RCSUtils::exportACARHSMovingAngles(solver, cfg, movingAngles);
		}

		const int timeCount = static_cast<int>(movingAngles.size());
		const int rowCount = solver.row;

		auto fillIncidentWave = [&](int timeIndex,
			double kinc[3], double einc[3], double hinc[3]) {
			const auto& angle = movingAngles[timeIndex];
			const double phiComputeDeg = isThetaPole(angle.thetaBodyDeg) ?
				canonicalPolePhiDeg() : angle.phiBodyDeg;
			EMSource pw;
			pw.initPW(data_pol, angle.thetaBodyDeg, phiComputeDeg);
			kinc[0] = pw.vW(0); kinc[1] = pw.vW(1); kinc[2] = pw.vW(2);
			einc[0] = pw.vE(0); einc[1] = pw.vE(1); einc[2] = pw.vE(2);
			hinc[0] = pw.vH(0); hinc[1] = pw.vH(1); hinc[2] = pw.vH(2);
			};

		auto computeColumn = [&](int timeIndex,
			std::vector<std::complex<double>>& rhs) {
			double kinc[3], einc[3], hinc[3];
			fillIncidentWave(timeIndex, kinc, einc, hinc);
			computeV(solver, kinc, einc, hinc);
			rhs.assign(solver.Vm.begin(), solver.Vm.end());
			};

		auto computeEntry = [&](int rowIndex, int timeIndex) {
			double kinc[3], einc[3], hinc[3];
			fillIncidentWave(timeIndex, kinc, einc, hinc);
			return computeVEntry(solver, kinc, einc, hinc, rowIndex);
			};

		ACARHS::Options options;
		options.tol = cfg.acaRHSTol;
		options.maxRank = cfg.acaRHSMaxRank;
		options.pivotEps = cfg.acaRHSPivotEps;

		std::cout << "ACA-RHS-MLFMA-Moving states = " << timeCount << std::endl;
		std::cout << "ACA-RHS-MLFMA-Moving RHS rows = " << rowCount << std::endl;
		auto acaRHSStart = std::chrono::high_resolution_clock::now();
		ACARHS::Factors factors =
			ACARHS::factorize(rowCount, timeCount, options, computeColumn, computeEntry);
		auto acaRHSEnd = std::chrono::high_resolution_clock::now();
		const auto acaRHSElapsed =
			std::chrono::duration_cast<std::chrono::microseconds>(acaRHSEnd - acaRHSStart);
		std::cout << "Info: ACA-RHS factorization done, elapsed: "
			<< acaRHSElapsed.count() / 1000.0 << " ms" << std::endl;
		std::cout << "ACA-RHS-MLFMA-Moving rank = " << factors.rank << std::endl;

		RCSUtils::exportACARHSCheck(cfg, movingAngles, factors, computeColumn);
		RCSUtils::exportACARHSFactors(cfg, factors);

		std::vector<std::complex<double>> modeSolutions(
			static_cast<size_t>(factors.rank) * rowCount,
			std::complex<double>(0.0, 0.0));

		for (int m = 0; m < factors.rank; ++m) {
			std::vector<std::complex<double>> rhsMode(rowCount);
			std::vector<std::complex<double>> xMode(
				rowCount, std::complex<double>(0.0, 0.0));

			for (int row = 0; row < rowCount; ++row) {
				rhsMode[row] = factors.U[static_cast<size_t>(m) * rowCount + row];
			}

			int itr_max = rowCount;
			int mr = 50;
			double tol_abs = 1.0e-8;
			double tol_rel = 1.0e-4;
			std::cout << "ACA-RHS-MLFMA-Moving solving mode "
				<< (m + 1) << " / " << factors.rank << std::endl;
			solver.matrix_solver(
				rowCount, xMode.data(), rhsMode.data(),
				itr_max, mr, tol_abs, tol_rel);

			for (int row = 0; row < rowCount; ++row) {
				modeSolutions[static_cast<size_t>(m) * rowCount + row] = xMode[row];
			}
		}

		auto rcsFile = cfg.openFile("_RCS_ACARHS_MLFMA_MOVING", "txt");
		if (!rcsFile.is_open()) {
			throw std::runtime_error("Cannot open ACA-RHS-MLFMA-Moving RCS file");
		}
		rcsFile << std::fixed << std::setprecision(12);
		rcsFile << "idx\ttime_s\trange_m\ttheta_global_deg\tphi_global_deg\t"
			<< "theta_body_deg\tphi_body_deg\taca_rank\trcs_db\trcs_v_db\trcs_h_db\n";

		std::vector<std::complex<double>> current(
			rowCount, std::complex<double>(0.0, 0.0));
		for (int timeIndex = 0; timeIndex < timeCount; ++timeIndex) {
			std::fill(current.begin(), current.end(), std::complex<double>(0.0, 0.0));
			for (int m = 0; m < factors.rank; ++m) {
				const std::complex<double> coeff =
					factors.V[static_cast<size_t>(m) * timeCount + timeIndex];
				for (int row = 0; row < rowCount; ++row) {
					current[row] += modeSolutions[static_cast<size_t>(m) * rowCount + row] * coeff;
				}
			}

			const auto& angle = movingAngles[timeIndex];
			const double phiComputeDeg = isThetaPole(angle.thetaBodyDeg) ?
				canonicalPolePhiDeg() : angle.phiBodyDeg;
			std::vector<std::complex<double>> RCS_pre(3, 0.0);
			RCSUtils::calculateRCS(
				solver, current.data(), RCS_pre,
				solver.wave.k1(), solver.wave.eta1(),
				angle.thetaBodyDeg, phiComputeDeg);
			RCSGridValue value = projectRCSGridValue(
				RCS_pre, data_pol, angle.thetaBodyDeg, phiComputeDeg);

			rcsFile << timeIndex << "\t" << angle.time_s << "\t" << angle.range << "\t"
				<< angle.thetaGlobalDeg << "\t" << angle.phiGlobalDeg << "\t"
				<< angle.thetaBodyDeg << "\t" << angle.phiBodyDeg << "\t"
				<< factors.rank << "\t"
				<< value.rcs << "\t" << value.rcs_v << "\t" << value.rcs_h << "\n";

			std::cout << "ACA-RHS-MLFMA-Moving time point "
				<< (timeIndex + 1) << " / " << timeCount << " RCS computed" << std::endl;
		}

		const size_t mem = solver.computeMem() +
			my_memory::memory1D<std::complex<double>>(factors.U) +
			my_memory::memory1D<std::complex<double>>(factors.V) +
			my_memory::memory1D<std::complex<double>>(modeSolutions);
		std::cout << "ACA-RHS-MLFMA-Moving rank = " << factors.rank << std::endl;
		std::cout << "Algorithm memory = "
			<< mem / (1024.0 * 1024.0) << " MB" << std::endl;
		std::cout << "RCS ACA-RHS-MLFMA-Moving export completed." << std::endl;
	}

	template<typename SolverType, typename ComputeVFunc>
	inline void computeMonoStatic_FSB_Moving(SolverType& solver, const RCSExportConfig& cfg,
		ComputeVFunc computeV, double rho = 1.2)
	{
		if (solver.integralEquType_ == 2) {
			throw std::runtime_error("FSB-MLFMA-Moving supports PEC mono only in this version.");
		}

		const std::string pol_wave = cfg.polWaveStr();
		const double data_pol = (pol_wave == "V") ? 0.0 : 90.0;
		const std::array<double, 3> radarPos = { cfg.radarX, cfg.radarY, cfg.radarZ };

		std::vector<MovingState> states = MovingTarget::readMovingPathFile(cfg.movingPathFile);
		std::vector<MovingAngleSample> movingAngles =
			MovingTarget::computeMovingAngleSamples(states, radarPos);
		if (cfg.exportMovingAngles) {
			RCSUtils::exportMovingAngles(solver, cfg, movingAngles);
		}

		std::vector<MovingAngleSample> uniqueAngles =
			MovingTarget::mergeMovingAngleSamples(movingAngles, cfg.movingAngleMergeTolDeg);

		const double radius = my_fsb::estimateBoundingSphereRadius(solver);
		const double W = std::max(1.0, rho * solver.wave.k1_abs() * radius);
		std::vector<std::vector<FSBMonoSample>> rings =
			my_fsb::buildMovingFSBSampleLayout(uniqueAngles, W, cfg.movingAngleMarginDeg);

		int expectedSampleCount = 0;
		for (const auto& ring : rings) {
			expectedSampleCount += static_cast<int>(ring.size());
		}
		if (expectedSampleCount <= 0) {
			throw std::runtime_error("FSB-MLFMA-Moving produced no FSB sample points.");
		}

		std::cout << "FSB-MLFMA-Moving states = " << states.size() << std::endl;
		std::cout << "FSB-MLFMA-Moving unique angle queries = " << uniqueAngles.size() << std::endl;
		std::cout << "FSB-MLFMA-Moving accurate sample points to compute = "
			<< expectedSampleCount << std::endl;

		int sampleCount = 0;
		for (auto& ring : rings) {
			for (auto& sample : ring) {
				EMSource pw;
				pw.initPW(data_pol, sample.thetaDeg, sample.phiDeg);
				double kinc[3] = { pw.vW(0), pw.vW(1), pw.vW(2) };
				double einc[3] = { pw.vE(0), pw.vE(1), pw.vE(2) };
				double hinc[3] = { pw.vH(0), pw.vH(1), pw.vH(2) };

				computeV(solver, kinc, einc, hinc);

				std::vector<std::complex<double>> I_solve(
					solver.row, std::complex<double>(0.0, 0.0));
				std::vector<std::complex<double>> Vec_R(solver.row);
				for (int i = 0; i < solver.row; ++i) {
					Vec_R[i] = solver.Vm[i];
				}

				int itr_max = solver.row;
				int mr = 50;
				double tol_abs = 1.0e-8;
				double tol_rel = 1.0e-4;
				solver.matrix_solver(
					solver.row, I_solve.data(), Vec_R.data(),
					itr_max, mr, tol_abs, tol_rel);

				std::vector<std::complex<double>> RCS_pre(3, 0.0);
				RCSUtils::calculateRCS(
					solver, I_solve.data(), RCS_pre,
					solver.wave.k1(), solver.wave.eta1(),
					sample.thetaDeg, sample.phiDeg);
				sample.field = std::move(RCS_pre);

				++sampleCount;
				std::cout << "FSB-MLFMA-Moving accurate sample point "
					<< sampleCount << " / " << expectedSampleCount
					<< " solved" << std::endl;
			}
		}

		if (cfg.exportMovingSamples) {
			RCSUtils::exportFSBSamples(solver, cfg, rings);
		}

		auto rcsFile = cfg.openFile("_RCS_FSB_MLFMA_MOVING", "txt");
		if (!rcsFile.is_open()) {
			throw std::runtime_error("Cannot open FSB-MLFMA-Moving RCS file");
		}
		std::ofstream fieldFile;
		if (cfg.exportMovingComplexField) {
			fieldFile = cfg.openFile("_FIELD_FSB_MLFMA_MOVING", "txt");
			if (!fieldFile.is_open()) {
				throw std::runtime_error("Cannot open FSB-MLFMA-Moving field file");
			}
			fieldFile << std::fixed << std::setprecision(12);
			fieldFile << "idx\ttime_s\trange_m\ttheta_body_deg\tphi_body_deg\t"
				<< "Ev_re\tEv_im\tEh_re\tEh_im\n";
		}

		rcsFile << std::fixed << std::setprecision(12);
		rcsFile << "idx\ttime_s\trange_m\ttheta_global_deg\tphi_global_deg\t"
			<< "theta_body_deg\tphi_body_deg\trcs_db\trcs_v_db\trcs_h_db\n";
		for (int i = 0; i < static_cast<int>(movingAngles.size()); ++i) {
			const auto& angle = movingAngles[i];
			const double phiInterpDeg = isThetaPole(angle.thetaBodyDeg) ?
				canonicalPolePhiDeg() : angle.phiBodyDeg;
			std::vector<std::complex<double>> RCS_pre = my_fsb::interpolateFSBField(
				rings, angle.thetaBodyDeg, phiInterpDeg);
			RCSGridValue value = projectRCSGridValue(
				RCS_pre, data_pol, angle.thetaBodyDeg, phiInterpDeg);

			rcsFile << i << "\t" << angle.time_s << "\t" << angle.range << "\t"
				<< angle.thetaGlobalDeg << "\t" << angle.phiGlobalDeg << "\t"
				<< angle.thetaBodyDeg << "\t" << angle.phiBodyDeg << "\t"
				<< value.rcs << "\t" << value.rcs_v << "\t" << value.rcs_h << "\n";

			if (cfg.exportMovingComplexField) {
				fieldFile << i << "\t" << angle.time_s << "\t" << angle.range << "\t"
					<< angle.thetaBodyDeg << "\t" << angle.phiBodyDeg << "\t"
					<< value.E_v.real() << "\t" << value.E_v.imag() << "\t"
					<< value.E_h.real() << "\t" << value.E_h.imag() << "\n";
			}
		}

		const size_t mem = solver.computeMem();
		std::cout << "FSB-MLFMA-Moving sample points = " << sampleCount << std::endl;
		std::cout << "Algorithm memory = "
			<< mem / (1024.0 * 1024.0) << " MB" << std::endl;
		std::cout << "RCS FSB-MLFMA-Moving export completed." << std::endl;
	}
	template<typename SolverType, typename ComputeVFunc>
	inline void computeMonoStatic_FSB(SolverType& solver, const RCSExportConfig& cfg,
		ComputeVFunc computeV, double rho = 1.2)
	{
		const std::string pol_wave = cfg.polWaveStr();
		const double data_pol = (pol_wave == "V") ? 0.0 : 90.0;
		const bool isPMCHWT = (solver.integralEquType_ == 2);
		const int totalRow = isPMCHWT ? 2 * solver.row : solver.row;
		const double radius = my_fsb::estimateBoundingSphereRadius(solver);
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
		std::cout << "FSB accurate sample points to compute = "
			<< expectedSampleCount << std::endl;

		std::vector<std::vector<FSBMonoSample>> rings;
		rings.reserve(thetaCount);

		int sampleCount = 0;
		for (int m = 1; m <= thetaCount; ++m) {
			const double theta = Pi * static_cast<double>(m) /
				static_cast<double>(thetaCount + 1);
			const double thetaDeg = theta * 180.0 / Pi;
			const int phiOrder = std::max(
				1, static_cast<int>(std::ceil(W * std::sin(theta))));
			const int phiCount = 2 * phiOrder + 1;

			std::vector<FSBMonoSample> ring;
			ring.reserve(phiCount);

			for (int n = -phiOrder; n <= phiOrder; ++n) {
				double phiDeg = my_fsb::normalizePhiDeg(
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
				std::cout << "FSB accurate sample point "
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

				FSBMonoSample sample;
				sample.thetaDeg = thetaDeg;
				sample.phiDeg = phiDeg;
				sample.phiOrder = phiOrder;
				sample.field = std::move(RCS_pre);
				ring.push_back(std::move(sample));
				++sampleCount;
			}

			rings.push_back(std::move(ring));
		}

		if (cfg.exportFSBSamples) {
			RCSUtils::exportFSBSamples(solver, cfg, rings);
		}
		RCSUtils::exportFSBScanGrid(solver, cfg, rings);

		const size_t mem = solver.computeMem();
		std::cout << "FSB sample points = " << sampleCount << std::endl;
		std::cout << "Algorithm memory = "
			<< mem / (1024.0 * 1024.0) << " MB" << std::endl;
		std::cout << "RCS FSB export completed." << std::endl;
	}
}
