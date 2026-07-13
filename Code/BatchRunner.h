#pragma once
#ifndef BATCHRUNNER_H
#define BATCHRUNNER_H

#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <iostream>
#include <stdexcept>
#include <cmath>
#include <exception>
#include <chrono>
#include <complex>
#include <omp.h>

#include "MoM.h"
#include "FMM.h"
#include "MLFMM.h"
#include "RWGGenerator.h"

struct MaterialParam {
	double epsilonR_real = 1.0;
	double epsilonR_imag = 0.0;
	double muR_real = 1.0;
	double muR_imag = 0.0;
};

struct SimTask {
	std::string nasFile;
	double freq;
	std::string polarization;
	double inc_th, inc_ph;
	double sca_th_s, sca_th_f;
	double sca_ph_s, sca_ph_f;
	double step;
	int selectAlgorithm;
	int selectIntegralEqu;
	int selectMatrixSolver;
	std::string selectMono_Dual;
	int N_points;
	double epsilonR_real, epsilonR_imag;
	double muR_real, muR_imag;
	int ompThreads;
};

struct BatchConfig {
	std::vector<std::string> nasFiles;
	std::vector<double> frequencies;
	std::vector<std::pair<double, double>> incidentAngles;

	double sca_th_s = 0.0, sca_th_f = 180.0;
	double sca_ph_s = 0.0, sca_ph_f = 360.0;
	double step = 1.0;

	std::vector<int> selectAlgorithms{ 2 };
	int selectIntegralEqu = 1;
	int selectMatrixSolver = 0;
	std::string selectMono_Dual = "dual";
	int N_points = 7;
	std::vector<MaterialParam> materialParams{ MaterialParam{} };
	int ompThreads = 12;

	std::vector<SimTask> generateTasks() const {
		std::vector<SimTask> tasks;

		for (const auto& nas : nasFiles) {
			for (double freq : frequencies) {
				for (const auto& [th, ph] : incidentAngles) {
					for (int algorithm : selectAlgorithms) {
						for (const auto& material : materialParams) {
							SimTask task;
							task.nasFile = nas;
							task.freq = freq;
							task.inc_th = th;
							task.inc_ph = ph;
							task.sca_th_s = sca_th_s;
							task.sca_th_f = sca_th_f;
							task.sca_ph_s = sca_ph_s;
							task.sca_ph_f = sca_ph_f;
							task.step = step;
							task.selectAlgorithm = algorithm;
							task.selectIntegralEqu = selectIntegralEqu;
							task.selectMatrixSolver = selectMatrixSolver;
							task.selectMono_Dual = selectMono_Dual;
							task.N_points = N_points;
							task.epsilonR_real = material.epsilonR_real;
							task.epsilonR_imag = material.epsilonR_imag;
							task.muR_real = material.muR_real;
							task.muR_imag = material.muR_imag;
							task.ompThreads = ompThreads;
							tasks.push_back(task);
						}
					}
				}
			}
		}

		return tasks;
	}
};

class BatchConfigParser {
public:
	static BatchConfig parse(const std::string& configPath) {
		BatchConfig cfg;

		std::ifstream file(configPath);
		if (!file.is_open()) {
			throw std::runtime_error("Cannot open config file: " + configPath);
		}

		bool hasGlobalMu = false;
		double globalMuReal = 1.0;
		double globalMuImag = 0.0;

		std::string line;
		while (std::getline(file, line)) {
			if (line.empty() || line[0] == '#') continue;

			auto pos = line.find('=');
			if (pos == std::string::npos) continue;

			std::string key = trim(line.substr(0, pos));
			std::string value = trim(line.substr(pos + 1));

			if (key == "nas_files") {
				cfg.nasFiles = splitString(value, ';');
			}
			else if (key == "frequencies") {
				cfg.frequencies.clear();
				for (const auto& s : splitString(value, ',')) {
					cfg.frequencies.push_back(std::stod(s));
				}
			}
			else if (key == "incident_angles") {
				cfg.incidentAngles.clear();
				for (const auto& pair : splitString(value, ',')) {
					auto parts = splitString(pair, ':');
					if (parts.size() == 2) {
						cfg.incidentAngles.emplace_back(std::stod(parts[0]), std::stod(parts[1]));
					}
				}
			}
			else if (key == "incident_theta_range") {
				auto parts = splitString(value, ':');
				if (parts.size() == 3) {
					double start = std::stod(parts[0]);
					double end = std::stod(parts[1]);
					double s = std::stod(parts[2]);
					int count = static_cast<int>(std::round((end - start) / s)) + 1;

					for (int i = 0; i < count; ++i) {
						double th = start + i * s;
						cfg.incidentAngles.emplace_back(th, 0.0);
					}
				}
			}
			else if (key == "incident_phi_range") {
				auto parts = splitString(value, ':');
				if (parts.size() == 3) {
					double start = std::stod(parts[0]);
					double end = std::stod(parts[1]);
					double s = std::stod(parts[2]);
					std::vector<std::pair<double, double>> expanded;

					if (cfg.incidentAngles.empty()) {
						for (double ph = start; ph <= end + 1e-9; ph += s) {
							expanded.emplace_back(90.0, ph);
						}
					}
					else {
						for (const auto& [th, oldPh] : cfg.incidentAngles) {
							(void)oldPh;
							for (double ph = start; ph <= end + 1e-9; ph += s) {
								expanded.emplace_back(th, ph);
							}
						}
					}

					cfg.incidentAngles = expanded;
				}
			}
			else if (key == "sca_theta") {
				auto parts = splitString(value, ':');
				if (parts.size() == 2) {
					cfg.sca_th_s = std::stod(parts[0]);
					cfg.sca_th_f = std::stod(parts[1]);
				}
			}
			else if (key == "sca_phi") {
				auto parts = splitString(value, ':');
				if (parts.size() == 2) {
					cfg.sca_ph_s = std::stod(parts[0]);
					cfg.sca_ph_f = std::stod(parts[1]);
				}
			}
			else if (key == "step") {
				cfg.step = std::stod(value);
			}
			else if (key == "algorithm") {
				cfg.selectAlgorithms.clear();
				for (const auto& s : splitString(value, ',')) {
					cfg.selectAlgorithms.push_back(std::stoi(s));
				}
				if (cfg.selectAlgorithms.empty()) {
					cfg.selectAlgorithms.push_back(2);
				}
			}
			else if (key == "integral_equation") {
				cfg.selectIntegralEqu = std::stoi(value);
			}
			else if (key == "matrix_solver") {
				cfg.selectMatrixSolver = std::stoi(value);
			}
			else if (key == "mono_dual") {
				cfg.selectMono_Dual = value;
			}
			else if (key == "gauss_points") {
				cfg.N_points = std::stoi(value);
			}
			else if (key == "epsilon_r") {
				cfg.materialParams.clear();

				MaterialParam material;
				auto parts = splitString(value, ',');
				if (parts.size() >= 1) material.epsilonR_real = std::stod(parts[0]);
				if (parts.size() >= 2) material.epsilonR_imag = std::stod(parts[1]);

				if (hasGlobalMu) {
					material.muR_real = globalMuReal;
					material.muR_imag = globalMuImag;
				}

				cfg.materialParams.push_back(material);
			}
			else if (key == "epsilon_r_list") {
				cfg.materialParams.clear();

				for (const auto& item : splitString(value, ';')) {
					auto parts = splitString(item, ',');
					if (parts.size() >= 1) {
						MaterialParam material;
						material.epsilonR_real = std::stod(parts[0]);
						if (parts.size() >= 2) material.epsilonR_imag = std::stod(parts[1]);

						if (hasGlobalMu) {
							material.muR_real = globalMuReal;
							material.muR_imag = globalMuImag;
						}

						cfg.materialParams.push_back(material);
					}
				}

				if (cfg.materialParams.empty()) {
					cfg.materialParams.push_back(MaterialParam{});
				}
			}
			else if (key == "mu_r") {
				auto parts = splitString(value, ',');
				if (parts.size() >= 1) globalMuReal = std::stod(parts[0]);
				if (parts.size() >= 2) globalMuImag = std::stod(parts[1]);
				hasGlobalMu = true;

				for (auto& material : cfg.materialParams) {
					material.muR_real = globalMuReal;
					material.muR_imag = globalMuImag;
				}
			}
			else if (key == "material_list") {
				cfg.materialParams.clear();

				for (const auto& item : splitString(value, ';')) {
					auto parts = splitString(item, ',');
					if (parts.size() >= 2) {
						MaterialParam material;
						material.epsilonR_real = std::stod(parts[0]);
						material.epsilonR_imag = std::stod(parts[1]);
						if (parts.size() >= 3) material.muR_real = std::stod(parts[2]);
						if (parts.size() >= 4) material.muR_imag = std::stod(parts[3]);

						if (hasGlobalMu && parts.size() < 3) {
							material.muR_real = globalMuReal;
							material.muR_imag = globalMuImag;
						}

						cfg.materialParams.push_back(material);
					}
				}

				if (cfg.materialParams.empty()) {
					cfg.materialParams.push_back(MaterialParam{});
				}
			}
			else if (key == "omp_threads") {
				cfg.ompThreads = std::stoi(value);
			}
		}

		return cfg;
	}

private:
	static std::string trim(const std::string& s) {
		size_t start = s.find_first_not_of(" \t\r\n");
		size_t end = s.find_last_not_of(" \t\r\n");
		return (start == std::string::npos) ? "" : s.substr(start, end - start + 1);
	}

	static std::vector<std::string> splitString(const std::string& s, char delim) {
		std::vector<std::string> result;
		std::stringstream ss(s);
		std::string item;

		while (std::getline(ss, item, delim)) {
			std::string trimmed = trim(item);
			if (!trimmed.empty()) {
				result.push_back(trimmed);
			}
		}

		return result;
	}
};

class BatchRunner {
public:
	static void run(const std::string& configPath) {
		std::cout << "========================================\n";
		std::cout << " Algorithm Batch Simulation Runner\n";
		std::cout << "========================================\n\n";

		BatchConfig batchCfg = BatchConfigParser::parse(configPath);
		std::vector<SimTask> tasks = batchCfg.generateTasks();

		std::cout << "Total tasks: " << tasks.size()
			<< " (x2 polarizations = " << tasks.size() * 2 << " runs)\n\n";

		auto batchStart = std::chrono::high_resolution_clock::now();

		for (size_t i = 0; i < tasks.size(); ++i) {
			const auto& task = tasks[i];

			std::cout << "----------------------------------------\n";
			std::cout << "[Task " << (i + 1) << "/" << tasks.size() << "]\n";
			std::cout << "  NAS: " << task.nasFile << "\n";
			std::cout << "  Freq: " << task.freq / 1e9 << " GHz\n";
			std::cout << "  Inc: theta=" << task.inc_th << ", phi=" << task.inc_ph << "\n";
			std::cout << "  Material: eps=(" << task.epsilonR_real << ", "
				<< task.epsilonR_imag << "), mu=("
				<< task.muR_real << ", " << task.muR_imag << ")\n";
			std::cout << "  Algo: " << task.selectAlgorithm
				<< ", IE: " << task.selectIntegralEqu
				<< ", Solver: " << task.selectMatrixSolver
				<< " (0=GMRES, 1=CGS)\n";
			std::cout << "  Polarization: h -> v (sequential)\n";
			std::cout << "----------------------------------------\n";

			auto taskStart = std::chrono::high_resolution_clock::now();
			int ret = executeSingleTask(task);
			auto taskEnd = std::chrono::high_resolution_clock::now();

			double elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(
				taskEnd - taskStart).count() / 1000.0;

			if (ret == 0) {
				std::cout << "[Task " << (i + 1) << "] DONE (" << elapsed << " s)\n\n";
			}
			else {
				std::cerr << "[Task " << (i + 1) << "] FAILED (code=" << ret << ")\n\n";
			}
		}

		auto batchEnd = std::chrono::high_resolution_clock::now();
		double totalTime = std::chrono::duration_cast<std::chrono::seconds>(
			batchEnd - batchStart).count();

		std::cout << "========================================\n";
		std::cout << "Batch complete. Total time: " << totalTime << " s\n";
		std::cout << "========================================\n";
	}

private:
	static int executeSingleTask(const SimTask& task) {
		try {
			omp_set_dynamic(0);
			omp_set_num_threads(task.ompThreads);

			if (task.selectMatrixSolver != 0 && task.selectMatrixSolver != 1) {
				throw std::invalid_argument(
					"Invalid matrix_solver value: "
					+ std::to_string(task.selectMatrixSolver)
					+ " (must be 0 for GMRES or 1 for CGS)");
			}

			std::vector<Point> points;
			std::vector<FaceElement> faces;

			if (readNasData(task.nasFile, points, faces) != 0) {
				std::cerr << "Error: Failed to read NAS file: " << task.nasFile << "\n";
				return 1;
			}

			std::vector<RWGBase> rwgBases;
			if (!RWGGenerator::GenerateRWG(faces, points, rwgBases)) {
				std::cerr << "Error: RWG generation failed\n";
				return 2;
			}

			EMSource waveForTree(task.freq, 90.0);
			std::complex<double> epsilonR(task.epsilonR_real, task.epsilonR_imag);
			std::complex<double> muR(task.muR_real, task.muR_imag);

			if (task.selectIntegralEqu == 2) {
				waveForTree.initDielectric(epsilonR, muR);
			}

			gaussPoints gausspoint(task.N_points);
			OCTree octree(rwgBases, waveForTree.wavelength(), task.nasFile);

			const std::string pols[2] = { "h", "v" };

			for (const auto& pol : pols) {
				std::cout << "  >> Computing polarization: " << pol << "\n";

				double polVal = (pol == "v") ? 0.0 : 90.0;
				EMSource wave(task.freq, polVal);

				if (task.selectIntegralEqu == 2) {
					wave.initDielectric(epsilonR, muR);
				}

				RCSExportConfig cfg(task.nasFile);

				if (task.selectAlgorithm < 0 || task.selectAlgorithm > 2) {
					throw std::invalid_argument(
						"Invalid algorithm value: "
						+ std::to_string(task.selectAlgorithm));
				}

				cfg.selectAlgo = static_cast<selectAlgorithm>(task.selectAlgorithm);
				cfg.material = (task.selectIntegralEqu == 2) ? materialType::Die : materialType::PEC;
				cfg.integralEqu = (task.selectIntegralEqu == 0) ? IntegralEquation::EFIE :
					(task.selectIntegralEqu == 2) ? IntegralEquation::PMCHWT :
					IntegralEquation::CFIE;
				cfg.monoDual = (task.selectMono_Dual == "mono") ? MonoDual::Mono : MonoDual::Dual;
				cfg.polWave = (pol == "v") ? Polarization::V : Polarization::H;
				cfg.nPoints = task.N_points;
				cfg.incTheta = task.inc_th;
				cfg.incPhi = task.inc_ph;
				cfg.scaThetaStart = task.sca_th_s;
				cfg.scaThetaEnd = task.sca_th_f;
				cfg.scaPhiStart = task.sca_ph_s;
				cfg.scaPhiEnd = task.sca_ph_f;
				cfg.scaStep = task.step;

				if (task.selectIntegralEqu == 2) {
					cfg.epsilonR = epsilonR;
					cfg.muR = muR;
				}

				cfg.validate();

				double E0 = 1.0;

				if (task.selectAlgorithm == 0) {
					MoM mom(cfg, task.selectIntegralEqu, task.selectMatrixSolver, task.selectMono_Dual,
						pol, octree.octreeNodes_, rwgBases,
						octree.maxLevel_, gausspoint, E0, wave);
				}
				else if (task.selectAlgorithm == 1) {
					FMM fmm(cfg, task.selectIntegralEqu, task.selectMatrixSolver, task.selectMono_Dual,
						pol, octree.octreeNodes_, rwgBases,
						octree.maxLevel_, gausspoint, E0, wave);
				}
				else if (task.selectAlgorithm == 2) {
					MLFMM mlfmm(cfg, task.selectIntegralEqu, task.selectMatrixSolver, task.selectMono_Dual,
						pol, octree.octreeNodes_, octree.octreeNodesDRvec_,
						rwgBases, octree.maxLevel_, gausspoint, E0, wave);
				}

				std::cout << "  >> Polarization " << pol << " done.\n";
			}

			return 0;
		}
		catch (const std::exception& e) {
			std::cerr << "Exception: " << e.what() << "\n";
			return -1;
		}
	}
};

#endif // BATCHRUNNER_H