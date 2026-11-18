#pragma once
#ifndef BATCHRUNNER_H
#define BATCHRUNNER_H

#include <algorithm>
#include <array>
#include <climits>
#include <cstddef>
#include <string>
#include <vector>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <limits>
#include <locale>
#include <sstream>
#include <iostream>
#include <stdexcept>
#include <cmath>
#include <exception>
#include <chrono>
#include <complex>
#include <optional>
#include <cctype>
#include <omp.h>

#include "OCTree.h"
#include "MoM.h"
#include "FMM.h"
#include "MLFMM.h"
#include "RWGGenerator.h"
#include "RunLogger.h"
#include "PathUtils.h"

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
	bool useFSBMono;
	bool useCalderonV1 = false;
	double fsbRho;
	bool exportFSBComplexField;
	bool exportFSBSamples;
	bool useFSBMLFMAMoving;
	std::string movingPathFile;
	double radarX, radarY, radarZ;
	double movingAngleMarginDeg;
	double movingAngleMergeTolDeg;
	bool exportMovingAngles;
	bool exportMovingComplexField;
	bool exportMovingSamples;
	bool useAccurateMLFMAMoving;
	bool useACARHSMoving;
	double acaRHSTol;
	int acaRHSMaxRank;
	double acaRHSPivotEps;
	bool exportACARHSCheck;
	int acaRHSCheckCount;
	bool exportACARHSFactors;
};

struct BatchConfig {
	std::vector<std::string> nasFiles;
	std::vector<double> frequencies;
	bool mergeRcsOutput = false;
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
	bool useFSBMono = false;
	bool useCalderonV1 = false;
	double fsbRho = 1.2;
	bool exportFSBComplexField = false;
	bool exportFSBSamples = false;
	bool useFSBMLFMAMoving = false;
	std::string movingPathFile;
	double radarX = 0.0, radarY = 0.0, radarZ = 0.0;
	double movingAngleMarginDeg = 2.0;
	double movingAngleMergeTolDeg = 0.25;
	bool exportMovingAngles = true;
	bool exportMovingComplexField = false;
	bool exportMovingSamples = false;
	bool useAccurateMLFMAMoving = false;
	bool useACARHSMoving = false;
	double acaRHSTol = 1.0e-3;
	int acaRHSMaxRank = 0;
	double acaRHSPivotEps = 1.0e-14;
	bool exportACARHSCheck = true;
	int acaRHSCheckCount = 10;
	bool exportACARHSFactors = false;

	std::vector<SimTask> generateTasks() const {
		std::vector<SimTask> tasks;
		const std::vector<MaterialParam> pecMaterial{ MaterialParam{} };
		const auto& activeMaterials = (selectIntegralEqu == 2) ? materialParams : pecMaterial;

		for (const auto& nas : nasFiles) {
			for (double freq : frequencies) {
				for (const auto& [th, ph] : incidentAngles) {
					for (int algorithm : selectAlgorithms) {
						for (const auto& material : activeMaterials) {
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
							task.useFSBMono = useFSBMono;
							task.useCalderonV1 = useCalderonV1;
							task.fsbRho = fsbRho;
							task.exportFSBComplexField = exportFSBComplexField;
							task.exportFSBSamples = exportFSBSamples;
							task.useFSBMLFMAMoving = useFSBMLFMAMoving;
							task.movingPathFile = movingPathFile;
							task.radarX = radarX;
							task.radarY = radarY;
							task.radarZ = radarZ;
							task.movingAngleMarginDeg = movingAngleMarginDeg;
							task.movingAngleMergeTolDeg = movingAngleMergeTolDeg;
							task.exportMovingAngles = exportMovingAngles;
							task.exportMovingComplexField = exportMovingComplexField;
							task.exportMovingSamples = exportMovingSamples;
							task.useAccurateMLFMAMoving = useAccurateMLFMAMoving;
							task.useACARHSMoving = useACARHSMoving;
							task.acaRHSTol = acaRHSTol;
							task.acaRHSMaxRank = acaRHSMaxRank;
							task.acaRHSPivotEps = acaRHSPivotEps;
							task.exportACARHSCheck = exportACARHSCheck;
							task.acaRHSCheckCount = acaRHSCheckCount;
							task.exportACARHSFactors = exportACARHSFactors;
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
		std::string frequencyKey;

		std::string line;
		std::size_t lineNo = 0;
		while (std::getline(file, line)) {
			++lineNo;
			if (line.empty() || line[0] == '#') continue;

			auto pos = line.find('=');
			if (pos == std::string::npos) continue;

			std::string key = trim(line.substr(0, pos));
			std::string value = trim(line.substr(pos + 1));

			if (key == "nas_files") {
				cfg.nasFiles = parseNasFiles(value, configPath, lineNo);
			}
			else if (key == "frequencies" || key == "frequency_range") {
				if (!frequencyKey.empty() && frequencyKey != key) {
					throw configValueError(key, value, lineNo,
						"frequencies and frequency_range are mutually exclusive");
				}
				frequencyKey = key;
				cfg.frequencies = key == "frequency_range"
					? parseFrequencyRange(value, key, lineNo)
					: parseFrequencyList(value, key, lineNo);
			}
			else if (key == "merge_rcs_output") {
				cfg.mergeRcsOutput = parseBool(value);
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
			else if (key == "incident_theta_range" || key == "incident_phi_range") {
				const auto [start, end, step] = parseRange(value, key, lineNo);
				const auto rangeError = [&](const char* message) {
					return configValueError(key, value, lineNo, message);
				};

				if (key == "incident_theta_range") {
					if (start < 0.0 || start > 180.0 || end < 0.0 || end > 180.0 ||
						step == 0.0 || (start < end && step < 0.0) || (start > end && step > 0.0)) {
						throw rangeError("expected theta in [0,180] and a nonzero step toward the end");
					}

					const double steps = std::round((end - start) / step);
					if (!std::isfinite(steps) || steps < 0.0 || steps > INT_MAX - 1) {
						throw rangeError("too many theta samples");
					}
					const double last = start + steps * step;
					if (last < 0.0 || last > 180.0) {
						throw rangeError("last theta sample is outside [0,180]");
					}

					const int count = static_cast<int>(steps) + 1;
					for (int i = 0; i < count; ++i) {
						cfg.incidentAngles.emplace_back(start + i * step, 0.0);
					}
				}
				else {
					if (start < 0.0 || end > 360.0 || start > end || step <= 0.0) {
						throw rangeError("expected ascending phi in [0,360] and a positive step");
					}

					const double limit = end + 1e-9;
					if (std::floor((limit - start) / step) > INT_MAX - 1) {
						throw rangeError("too many phi samples");
					}

					std::vector<double> phiValues;
					for (double phi = start; phi <= limit;) {
						if (phiValues.size() == static_cast<std::size_t>(INT_MAX)) {
							throw rangeError("too many phi samples");
						}
						phiValues.push_back(phi);
						const double next = phi + step;
						if (next <= phi) {
							throw rangeError("phi step does not advance the angle");
						}
						phi = next;
					}

					if (cfg.incidentAngles.empty()) {
						cfg.incidentAngles.emplace_back(90.0, 0.0);
					}
					std::vector<std::pair<double, double>> expanded;
					if (cfg.incidentAngles.size() > expanded.max_size() / phiValues.size()) {
						throw rangeError("too many incident angles");
					}
					expanded.reserve(cfg.incidentAngles.size() * phiValues.size());
					for (const auto& angle : cfg.incidentAngles) {
						for (double phi : phiValues) {
							expanded.emplace_back(angle.first, phi);
						}
					}
					cfg.incidentAngles.swap(expanded);
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
					cfg.selectAlgorithms.push_back(parseSelection(s, key));
				}
				if (cfg.selectAlgorithms.empty()) {
					cfg.selectAlgorithms.push_back(2);
				}
			}
			else if (key == "integral_equation") {
				cfg.selectIntegralEqu = parseSelection(value, key);
			}
			else if (key == "matrix_solver") {
				cfg.selectMatrixSolver = parseSelection(value, key);
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
			else if (key == "use_fsb_mono") {
				cfg.useFSBMono = parseBool(value);
			}
			else if (key == "calderon_v1") {
				cfg.useCalderonV1 = parseBool(value);
			}
			else if (key == "fsb_rho") {
				cfg.fsbRho = std::stod(value);
			}
			else if (key == "export_fsb_complex_field") {
				cfg.exportFSBComplexField = parseBool(value);
			}
			else if (key == "export_fsb_samples") {
				cfg.exportFSBSamples = parseBool(value);
			}
			else if (key == "use_fsb_mlfma_moving") {
				cfg.useFSBMLFMAMoving = parseBool(value);
			}
			else if (key == "moving_path_file") {
				cfg.movingPathFile = value;
			}
			else if (key == "radar_position") {
				auto parts = splitString(value, ',');
				if (parts.size() == 3) {
					cfg.radarX = std::stod(parts[0]);
					cfg.radarY = std::stod(parts[1]);
					cfg.radarZ = std::stod(parts[2]);
				}
			}
			else if (key == "moving_angle_margin_deg") {
				cfg.movingAngleMarginDeg = std::stod(value);
			}
			else if (key == "moving_angle_merge_tol_deg") {
				cfg.movingAngleMergeTolDeg = std::stod(value);
			}
			else if (key == "export_moving_angles") {
				cfg.exportMovingAngles = parseBool(value);
			}
			else if (key == "export_moving_complex_field") {
				cfg.exportMovingComplexField = parseBool(value);
			}
			else if (key == "export_moving_samples") {
				cfg.exportMovingSamples = parseBool(value);
			}
			else if (key == "use_accurate_mlfma_moving") {
				cfg.useAccurateMLFMAMoving = parseBool(value);
			}
			else if (key == "use_aca_rhs_moving") {
				cfg.useACARHSMoving = parseBool(value);
			}
			else if (key == "aca_rhs_tol") {
				cfg.acaRHSTol = std::stod(value);
			}
			else if (key == "aca_rhs_max_rank") {
				cfg.acaRHSMaxRank = std::stoi(value);
			}
			else if (key == "aca_rhs_pivot_eps") {
				cfg.acaRHSPivotEps = std::stod(value);
			}
			else if (key == "export_aca_rhs_check") {
				cfg.exportACARHSCheck = parseBool(value);
			}
			else if (key == "aca_rhs_check_count") {
				cfg.acaRHSCheckCount = std::stoi(value);
			}
			else if (key == "export_aca_rhs_factors") {
				cfg.exportACARHSFactors = parseBool(value);
			}
		}

		return cfg;
	}

private:
	static std::invalid_argument configValueError(const std::string& key,
		const std::string& value, std::size_t lineNo, const char* message) {
		return std::invalid_argument("Line " + std::to_string(lineNo) + ": " + key + " = " + value + " (" + message + ")");
	}

	static std::vector<std::string> parseNasFiles(const std::string& value,
		const std::string& configPath, std::size_t lineNo) {
		namespace fs = std::filesystem;
		auto files = splitString(value, ';');
		fs::path configDirectory;
		for (auto& entry : files) {
			const auto pathError = [&](const std::string& message) {
				return configValueError("nas_files", entry, lineNo, message.c_str());
			};
			try {
				const fs::path input(entry);
				if (input.is_absolute()) continue;

				const fs::path parent = input.parent_path();
				if (input.has_root_path() || (!parent.empty() && parent != ".") ||
					entry.find_first_of(":*?") != std::string::npos) {
					throw pathError("expected a NAS filename in the batch configuration directory");
				}
				std::string extension = input.extension().string();
				for (char& ch : extension) {
					ch = static_cast<char>(std::tolower(static_cast<unsigned char>(ch)));
				}
				if (extension != ".nas") throw pathError("expected a .nas file");

				if (configDirectory.empty()) {
					configDirectory = fs::canonical(fs::absolute(configPath).parent_path());
				}
				const fs::path candidate = configDirectory / input.filename();
				const fs::path resolved = fs::canonical(candidate);
				if (!fs::is_regular_file(resolved)) {
					throw pathError("expected a regular NAS file at " + candidate.string());
				}
				if (!fs::equivalent(resolved.parent_path(), configDirectory)) {
					throw pathError("NAS path resolves outside the batch configuration directory: " +
						candidate.string() + " -> " + resolved.string());
				}
				entry = resolved.string();
			}
			catch (const fs::filesystem_error& e) {
				throw pathError(e.what());
			}
		}
		return files;
	}

	static bool parseFiniteNumber(const std::string& field, double& number) {
		std::size_t consumed = 0;
		try {
			number = std::stod(field, &consumed);
		}
		catch (const std::logic_error&) {
			return false;
		}
		return consumed == field.size() && std::isfinite(number);
	}

	static std::array<double, 3> parseRange(const std::string& value,
		const std::string& key, std::size_t lineNo) {
		const auto first = value.find(':');
		const auto second = first == std::string::npos ? first : value.find(':', first + 1);
		if (second == std::string::npos || value.find(':', second + 1) != std::string::npos) {
			throw configValueError(key, value, lineNo, "expected start:end:step");
		}

		const std::array<std::string, 3> fields{
			value.substr(0, first),
			value.substr(first + 1, second - first - 1),
			value.substr(second + 1)
		};
		std::array<double, 3> range{};
		for (std::size_t i = 0; i < fields.size(); ++i) {
			if (!parseFiniteNumber(trim(fields[i]), range[i])) {
				throw configValueError(key, value, lineNo, "expected three finite numbers");
			}
		}
		return range;
	}

	static std::vector<double> parseFrequencyList(const std::string& value,
		const std::string& key, std::size_t lineNo) {
		const auto fields = splitString(value, ',');
		if (fields.empty()) {
			throw configValueError(key, value, lineNo, "expected at least one frequency in Hz");
		}
		std::vector<double> frequencies;
		frequencies.reserve(fields.size());
		for (const auto& field : fields) {
			double frequency = 0.0;
			if (!parseFiniteNumber(field, frequency) || frequency <= 0.0) {
				throw configValueError(key, value, lineNo, "expected positive finite frequencies in Hz");
			}
			frequencies.push_back(frequency);
		}
		return frequencies;
	}

	static std::vector<double> parseFrequencyRange(const std::string& value,
		const std::string& key, std::size_t lineNo) {
		const auto [start, end, step] = parseRange(value, key, lineNo);
		if (start <= 0.0 || end < start || step <= 0.0) {
			throw configValueError(key, value, lineNo, "expected 0 < start <= end and step > 0 in Hz");
		}

		constexpr std::size_t maxSweepPoints = 100'000;
		std::vector<double> frequencies{start};
		if (start == end) return frequencies;

		const double endSpacing = end - std::nextafter(end, 0.0);
		const double tolerance = std::min(2.0 * endSpacing, step / 8.0);
		for (std::size_t i = 1; ; ++i) {
			const double frequency = std::fma(static_cast<double>(i), step, start);
			const double delta = frequency - end;
			if (delta > tolerance) break;
			if (frequency <= frequencies.back()) {
				throw configValueError(key, value, lineNo, "frequency step does not advance the frequency");
			}
			if (frequencies.size() == maxSweepPoints) {
				throw configValueError(key, value, lineNo, "too many frequency points (maximum 100000)");
			}
			if (std::abs(delta) <= tolerance) {
				frequencies.push_back(end);
				break;
			}
			frequencies.push_back(frequency);
		}
		return frequencies;
	}

	static int parseSelection(const std::string& value, const std::string& key) {
		std::size_t consumed = 0;
		try {
			const int selection = std::stoi(value, &consumed);
			if (consumed == value.size()) return selection;
		}
		catch (const std::logic_error&) {
			throw std::invalid_argument("Invalid " + key + " value: " + value);
		}
		throw std::invalid_argument("Invalid " + key + " value: " + value);
	}

	static bool parseBool(const std::string& s) {
		std::string v = trim(s);
		for (char& ch : v) {
			ch = static_cast<char>(std::tolower(static_cast<unsigned char>(ch)));
		}
		return v == "1" || v == "true" || v == "yes" || v == "on";
	}

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
		BatchConfig batchCfg = BatchConfigParser::parse(configPath);
		std::vector<SimTask> tasks = batchCfg.generateTasks();
		std::optional<ScopedRunLogger> runLogger;
		std::ofstream sweepOutput;
		std::filesystem::path sweepPath, partialPath;

		std::cout << "========================================\n";
		std::cout << " Algorithm Batch Simulation Runner\n";
		std::cout << "========================================\n\n";
		if (batchCfg.mergeRcsOutput) {
			if (batchCfg.selectMono_Dual != "mono" || batchCfg.selectIntegralEqu != 1 ||
				batchCfg.useFSBMono || batchCfg.useFSBMLFMAMoving ||
				batchCfg.useAccurateMLFMAMoving || batchCfg.useACARHSMoving || batchCfg.useCalderonV1 ||
				std::any_of(batchCfg.selectAlgorithms.begin(), batchCfg.selectAlgorithms.end(),
					[](int algorithm) { return algorithm != 2; })) {
				throw std::invalid_argument("merge_rcs_output requires ordinary static monostatic MLFMA-CFIE");
			}
			if (tasks.empty()) throw std::invalid_argument("merge_rcs_output requires at least one simulation task");

			const std::filesystem::path directory(getParentPath(configPath));
			const std::string stem = getBaseName(configPath) + "_RCS_SWEEP_" + getCurrentTimeString();
			std::string name = stem;
			for (std::size_t suffix = 2; ; ++suffix) {
				sweepPath = directory / (name + ".txt");
				partialPath = directory / (name + ".part");
				if (!std::filesystem::exists(sweepPath) && !std::filesystem::exists(partialPath)) break;
				name = stem + "_" + std::to_string(suffix);
			}
			sweepOutput.open(partialPath);
			if (!sweepOutput) throw std::runtime_error("Failed to open RCS sweep output: " + partialPath.string());
			sweepOutput.imbue(std::locale::classic());
			sweepOutput << "task_id\tfrequency_hz\tincident_pol\ttheta_deg\tphi_deg"
				"\trcs_total_dbsm\trcs_v_dbsm\trcs_h_dbsm\n";
			sweepOutput.flush();
			if (!sweepOutput) throw std::runtime_error("Failed to write RCS sweep header: " + partialPath.string());
			std::cout << "RCS sweep output: " << sweepPath.string() << "\n\n";
		}
		else {
			runLogger.emplace(makeBatchLogPath(configPath));
			if (runLogger->enabled()) {
				std::cout << "Run log file: " << runLogger->path() << "\n\n";
			}
		}

		auto* mergedRcs = batchCfg.mergeRcsOutput ? &sweepOutput : nullptr;
		const bool numberOutputTasks = batchCfg.mergeRcsOutput || batchCfg.frequencies.size() > 1;

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
			std::cout << "  FSB Mono: " << (task.useFSBMono ? "on" : "off")
				<< ", rho=" << task.fsbRho
				<< ", field=" << (task.exportFSBComplexField ? "on" : "off")
				<< ", samples=" << (task.exportFSBSamples ? "on" : "off")
				<< ", moving=" << (task.useFSBMLFMAMoving ? "on" : "off") << "\n";
			std::cout << "  Accurate MLFMA Moving: "
				<< (task.useAccurateMLFMAMoving ? "on" : "off") << "\n";
			std::cout << "  ACA-RHS Moving: " << (task.useACARHSMoving ? "on" : "off")
				<< ", tol=" << task.acaRHSTol
				<< ", maxRank=" << task.acaRHSMaxRank << "\n";
			std::cout << "  Polarization: h -> v (sequential)\n";
			std::cout << "----------------------------------------\n";

			auto taskStart = std::chrono::high_resolution_clock::now();
			int ret = executeSingleTask(task, numberOutputTasks ? i + 1 : 0, mergedRcs);
			if (mergedRcs) {
				mergedRcs->flush();
				if (ret != 0 || !*mergedRcs) {
					std::ostringstream error;
					error << "RCS sweep failed at task " << (i + 1) << " ("
						<< std::setprecision(std::numeric_limits<double>::max_digits10)
						<< task.freq << " Hz";
					if (ret != 0) error << ", code " << ret;
					else error << ", write error";
					error << "). Partial results: " << partialPath.string();
					throw std::runtime_error(error.str());
				}
			}
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

		if (mergedRcs) {
			sweepOutput.close();
			if (!sweepOutput) throw std::runtime_error("Failed to close RCS sweep output: " + partialPath.string());
			std::filesystem::rename(partialPath, sweepPath);
			std::cout << "RCS sweep saved: " << sweepPath.string() << "\n";
		}

		auto batchEnd = std::chrono::high_resolution_clock::now();
		double totalTime = std::chrono::duration_cast<std::chrono::seconds>(
			batchEnd - batchStart).count();

		std::cout << "========================================\n";
		std::cout << "Batch complete. Total time: " << totalTime << " s\n";
		std::cout << "========================================\n";
	}

private:
	static std::string makeBatchLogPath(const std::string& configPath) {
		return getParentPath(configPath) + "\\" + getBaseName(configPath) +
			"_RUN_LOG_" + getCurrentTimeString() + ".txt";
	}

	static int executeSingleTask(const SimTask& task, std::size_t outputTaskIndex, std::ostream* mergedRcs) {
		try {
			omp_set_dynamic(0);
			omp_set_num_threads(task.ompThreads);

			if (task.selectMono_Dual != "mono" && task.selectMono_Dual != "dual") {
				throw std::invalid_argument("Invalid mono_dual value: " + task.selectMono_Dual);
			}

			RCSExportConfig cfg(task.nasFile);
			if (mergedRcs) {
				cfg.mergedRcs = mergedRcs;
				cfg.sweepTaskIndex = outputTaskIndex;
				cfg.sweepFrequencyHz = task.freq;
			}
			else if (outputTaskIndex != 0) {
				std::ostringstream tag;
				tag.imbue(std::locale::classic());
				tag << "_f" << std::defaultfloat << std::setprecision(std::numeric_limits<double>::max_digits10)
					<< task.freq << "Hz_t" << outputTaskIndex;
				cfg.batchOutputTag = tag.str();
			}
			cfg.selectAlgo = static_cast<selectAlgorithm>(task.selectAlgorithm);
			cfg.integralEqu = static_cast<IntegralEquation>(task.selectIntegralEqu);
			cfg.monoDual = (task.selectMono_Dual == "mono") ? MonoDual::Mono : MonoDual::Dual;
			cfg.nPoints = task.N_points;
			cfg.incTheta = task.inc_th;
			cfg.incPhi = task.inc_ph;
			cfg.scaThetaStart = task.sca_th_s;
			cfg.scaThetaEnd = task.sca_th_f;
			cfg.scaPhiStart = task.sca_ph_s;
			cfg.scaPhiEnd = task.sca_ph_f;
			cfg.scaStep = task.step;
			cfg.useFSBMono = task.useFSBMono;
			cfg.useCalderonV1 = task.useCalderonV1;
			cfg.fsbRho = task.fsbRho;
			cfg.exportFSBComplexField = task.exportFSBComplexField;
			cfg.exportFSBSamples = task.exportFSBSamples;
			cfg.useFSBMLFMAMoving = task.useFSBMLFMAMoving;
			cfg.movingPathFile = task.movingPathFile;
			cfg.radarX = task.radarX;
			cfg.radarY = task.radarY;
			cfg.radarZ = task.radarZ;
			cfg.movingAngleMarginDeg = task.movingAngleMarginDeg;
			cfg.movingAngleMergeTolDeg = task.movingAngleMergeTolDeg;
			cfg.exportMovingAngles = task.exportMovingAngles;
			cfg.exportMovingComplexField = task.exportMovingComplexField;
			cfg.exportMovingSamples = task.exportMovingSamples;
			cfg.useAccurateMLFMAMoving = task.useAccurateMLFMAMoving;
			cfg.useACARHSMoving = task.useACARHSMoving;
			cfg.acaRHSTol = task.acaRHSTol;
			cfg.acaRHSMaxRank = task.acaRHSMaxRank;
			cfg.acaRHSPivotEps = task.acaRHSPivotEps;
			cfg.exportACARHSCheck = task.exportACARHSCheck;
			cfg.acaRHSCheckCount = task.acaRHSCheckCount;
			cfg.exportACARHSFactors = task.exportACARHSFactors;
			if (cfg.integralEqu == IntegralEquation::PMCHWT) {
				cfg.epsilonR = {task.epsilonR_real, task.epsilonR_imag};
				cfg.muR = {task.muR_real, task.muR_imag};
			}
			cfg.finalizeTask(task.selectMatrixSolver, task.freq);
			gaussPoints gausspoint(cfg.nPoints);

			std::vector<Point> points;
			std::vector<FaceElement> faces;

			MeshReadReport meshReport;
			if (readNasData(task.nasFile, points, faces, meshReport) != 0) {
				std::cerr << "Error: Failed to read NAS file: " << task.nasFile << "\n";
				return 1;
			}

			std::vector<RWGBase> rwgBases;
			RWGOptions rwgOptions;
			rwgOptions.requireClosedSurface = cfg.integralEqu != IntegralEquation::EFIE;
			RWGGenerationReport rwgReport;
			if (!RWGGenerator::GenerateRWG(faces, points, rwgBases, rwgOptions, rwgReport)) {
				std::cerr << "Error: RWG generation failed\n";
				return 2;
			}

			std::optional<CalderonGeometry> cpGeometry;
			const CalderonGeometry* cpData = nullptr;
			if (cfg.useCalderonV1)
				cpData = &cpGeometry.emplace(points, faces, rwgBases, rwgReport);

			std::optional<OCTree> octree;
			if (cfg.selectAlgo != selectAlgorithm::MoM) {
				EMSource waveForTree(task.freq, 90.0);
				if (cfg.material == materialType::Die) {
					waveForTree.initDielectric(cfg.epsilonR, cfg.muR);
				}
				octree.emplace(rwgBases, waveForTree.wavelength(), task.nasFile);
			}

			const double E0 = 1.0;
			for (const auto pol : {Polarization::H, Polarization::V}) {
				cfg.polWave = pol;
				const bool vertical = pol == Polarization::V;
				const char* polName = vertical ? "v" : "h";
				std::cout << "  >> Computing polarization: " << polName << "\n";
				EMSource wave(task.freq, vertical ? 0.0 : 90.0);
				if (cfg.material == materialType::Die) {
					wave.initDielectric(cfg.epsilonR, cfg.muR);
				}

				if (cfg.selectAlgo == selectAlgorithm::MoM) {
					MoM mom(cfg, task.selectMatrixSolver, rwgBases, gausspoint, E0, wave);
				}
				else if (cfg.selectAlgo == selectAlgorithm::FMM) {
					FMM fmm(cfg, task.selectMatrixSolver, octree->octreeNodes_, rwgBases,
						octree->maxLevel_, gausspoint, E0, wave);
				}
				else if (cfg.selectAlgo == selectAlgorithm::MLFMA) {
					MLFMM mlfmm(cfg, task.selectMatrixSolver, octree->octreeNodes_, octree->octreeNodesDRvec_,
						rwgBases, octree->maxLevel_, gausspoint, E0, wave, cpData);
				}

				std::cout << "  >> Polarization " << polName << " done.\n";
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
