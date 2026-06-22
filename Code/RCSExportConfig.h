// RCSExportConfig.h
#pragma once
#ifndef RCS_EXPORT_CONFIG_H
#define RCS_EXPORT_CONFIG_H
#include <complex>
#include <filesystem>
#include <fstream>
#include <string>
#include <stdexcept>
#include <sstream>
#include <iomanip>

#include "PathUtils.h"

enum class selectAlgorithm { MoM = 0, FMM = 1, MLFMA = 2 };
enum class materialType { PEC = 0, Die = 1 };
enum class IntegralEquation { EFIE = 0, CFIE = 1, PMCHWT = 2 };
enum class MonoDual { Mono, Dual };
enum class Polarization { V, H };

struct RCSExportConfig {
	selectAlgorithm selectAlgo{ selectAlgorithm::MoM };
	materialType material{ materialType::PEC };
	IntegralEquation integralEqu{ IntegralEquation::EFIE };
	MonoDual monoDual{ MonoDual::Dual };
	Polarization polWave{ Polarization::V };
	int nPoints{ 0 };
	double incTheta{ 0.0 }, incPhi{ 0.0 };
	double scaThetaStart{ 0.0 }, scaPhiStart{ 0.0 };
	double scaThetaEnd{ 0.0 }, scaPhiEnd{ 0.0 };
	double scaStep{ 1.0 };
	std::complex<double> epsilonR{ 1.0,0.0 };
	std::complex<double> muR{ 1.0,0.0 };

	std::string nasPath;
	std::string parentPath;
	std::string baseName;
	std::string fixedTimeStamp;

	// 默认构造函数
	RCSExportConfig() = default;

	explicit RCSExportConfig(const std::string& path) : nasPath(path) {
		updatePathInfo();
	}

	// 完整构造函数
	RCSExportConfig(
		const std::string& path,
		IntegralEquation ie, MonoDual md, Polarization pw, int np,
		double ith, double iph, double sts, double sps,
		double ste, double spe, double step = 1.0)
		: nasPath(path), integralEqu(ie), monoDual(md), polWave(pw), nPoints(np),
		incTheta(ith), incPhi(iph), scaThetaStart(sts), scaPhiStart(sps),
		scaThetaEnd(ste), scaPhiEnd(spe), scaStep(step) {
		updatePathInfo();
		validate();
	}

	// 当 nasPath 改变时，重新计算 parentPath 和 baseName
	void setNasPath(const std::string& path) {
		nasPath = path;
		updatePathInfo();
	}

	// 获取时间戳
	void setFixedTimeStamp(const std::string& timestamp) {
		fixedTimeStamp = timestamp;
	}

	std::string getTimeStamp() const {
		return fixedTimeStamp.empty() ? getCurrentTimeString() : fixedTimeStamp;
	}

	// 验证参数是否有效
	bool isValid() const noexcept {
		if (nPoints <= 0) return false;
		if (scaStep <= 0.0) return false;
		if (incTheta < 0.0 || incTheta > 180.0) return false;
		if (incPhi < 0.0 || incPhi > 360.0) return false;
		if (scaThetaStart < 0.0 || scaThetaStart > 180.0) return false;
		if (scaThetaEnd < 0.0 || scaThetaEnd > 180.0) return false;
		if (scaPhiStart < 0.0 || scaPhiStart > 360.0) return false;
		if (scaPhiEnd < 0.0 || scaPhiEnd > 360.0) return false;
		return true;
	}

	void validate() const {
		if (!isValid()) {
			throw std::invalid_argument("RCSExportConfig: invalid parameters");
		}
	}

	// 枚举转字符串（用于文件名）
	std::string algoStr() const {
		switch (selectAlgo) {
		case selectAlgorithm::MoM:   return "MoM";
		case selectAlgorithm::FMM:   return "FMM";
		case selectAlgorithm::MLFMA: return "MLFMA";
		default: return "Unknown";
		}
	}

	std::string materialStr() const {
		switch (material) {
		case materialType::PEC:   return "PEC";
		case materialType::Die:   return "Die";
		default: return "Unknown";
		}
	}

	std::string integralEquStr() const {
		switch (integralEqu) {
		case IntegralEquation::EFIE: return "EFIE";
		case IntegralEquation::CFIE: return "CFIE";
		case IntegralEquation::PMCHWT: return "PMCHWT";

		default: return "Unknown";
		}
	}

	std::string monoDualStr() const {
		return (monoDual == MonoDual::Mono) ? "Mono" : "Dual";
	}

	std::string polWaveStr() const {
		return (polWave == Polarization::V) ? "V" : "H";
	}

	// 生成文件名的核心部分（不含路径、前缀、时间戳）
	std::string toCoreString() const {
		std::ostringstream oss;
		oss << algoStr() << "_"
			<< materialStr() << "_"
			<< integralEquStr() << "_"
			<< monoDualStr() << "_"
			<< polWaveStr() << "_"
			<< nPoints << "pts";
		if (material == materialType::Die) {
			oss << "_eps" << epsilonR.real() << "_" << epsilonR.imag();
		}
		if (monoDual == MonoDual::Dual) {
			oss << "_thetai" << static_cast<int>(incTheta)
				<< "_phii" << static_cast<int>(incPhi);
		}
		oss << "_thetas" << static_cast<int>(scaThetaStart) << "_"
			<< static_cast<int>(scaThetaEnd)
			<< "_phis" << static_cast<int>(scaPhiStart) << "_"
			<< static_cast<int>(scaPhiEnd)
			<< "_step" << std::fixed << std::setprecision(1) << scaStep;
		return oss.str();
	}
	// ----- 通用文件名生成和文件打开 -----
	// 构建完整路径：基础名 + 核心串 + 自定义后缀 + 时间戳 + 扩展名
	std::string buildFilePath(const std::string& suffix, const std::string& extension = "txt") const {
		namespace fs = std::filesystem;
		std::string time = getTimeStamp();
		std::string filename = baseName + "_" + toCoreString() + suffix + "_" + time + "." + extension;
		return (fs::path(parentPath) / filename).string();
	}
	// 直接打开文件，返回 ofstream (mode 默认参数为std::ios::out => 文本模式|  std::ios::binary => 二进制模式)
	std::ofstream openFile(
		const std::string& suffix,
		const std::string& extension = "txt",
		std::ios::openmode mode = std::ios::out) const {
		namespace fs = std::filesystem;
		fs::create_directories(parentPath);
		std::string fullPath = buildFilePath(suffix, extension);
		std::ofstream file(fullPath, mode);
		if (!file.is_open()) {
			throw std::runtime_error("Failed to open file: " + fullPath);
		}
		return file;
	}
private:
	void updatePathInfo() {
		parentPath = getParentPath(nasPath);
		baseName = getBaseName(nasPath);
	}
};
#endif // RCS_EXPORT_CONFIG_H