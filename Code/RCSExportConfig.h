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
#include <cmath>

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
	bool useFSBMono{ false };
	double fsbRho{ 1.2 };
	bool exportFSBComplexField{ false };
	bool exportFSBSamples{ false };
	bool useFSBMLFMAMoving{ false };
	std::string movingPathFile;
	double radarX{ 0.0 };
	double radarY{ 0.0 };
	double radarZ{ 0.0 };
	double movingAngleMarginDeg{ 2.0 };
	double movingAngleMergeTolDeg{ 0.25 };
	bool exportMovingAngles{ true };
	bool exportMovingComplexField{ false };
	bool exportMovingSamples{ false };
	bool useAccurateMLFMAMoving{ false };
	bool useACARHSMoving{ false };
	double acaRHSTol{ 1.0e-3 };
	int acaRHSMaxRank{ 0 };
	double acaRHSPivotEps{ 1.0e-14 };
	bool exportACARHSCheck{ true };
	int acaRHSCheckCount{ 10 };
	bool exportACARHSFactors{ false };
	std::complex<double> epsilonR{ 1.0,0.0 };
	std::complex<double> muR{ 1.0,0.0 };

	std::string nasPath;
	std::string parentPath;
	std::string baseName;
	std::string fixedTimeStamp;

	// Ĭ�Ϲ��캯��
	RCSExportConfig() = default;

	explicit RCSExportConfig(const std::string& path) : nasPath(path) {
		updatePathInfo();
	}

	// �������캯��
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

	// �� nasPath �ı�ʱ�����¼��� parentPath �� baseName
	void setNasPath(const std::string& path) {
		nasPath = path;
		updatePathInfo();
	}

	// ��ȡʱ���
	void setFixedTimeStamp(const std::string& timestamp) {
		fixedTimeStamp = timestamp;
	}

	std::string getTimeStamp() const {
		return fixedTimeStamp.empty() ? getCurrentTimeString() : fixedTimeStamp;
	}

	// ��֤�����Ƿ���Ч
	bool isValid() const noexcept {
		if (nPoints <= 0) return false;
		if (scaStep <= 0.0) return false;
		if (incTheta < 0.0 || incTheta > 180.0) return false;
		if (incPhi < 0.0 || incPhi > 360.0) return false;
		if (scaThetaStart < 0.0 || scaThetaStart > 180.0) return false;
		if (scaThetaEnd < 0.0 || scaThetaEnd > 180.0) return false;
		if (scaPhiStart < 0.0 || scaPhiStart > 360.0) return false;
		if (scaPhiEnd < 0.0 || scaPhiEnd > 360.0) return false;
		if (fsbRho <= 0.0) return false;
		const int movingModeCount =
			(useAccurateMLFMAMoving ? 1 : 0) +
			(useACARHSMoving ? 1 : 0) +
			(useFSBMLFMAMoving ? 1 : 0);
		if (movingModeCount > 1) return false;
		const bool hasMovingRadar =
			(radarX * radarX + radarY * radarY + radarZ * radarZ) > 0.0;
		if (useAccurateMLFMAMoving) {
			if (selectAlgo != selectAlgorithm::MLFMA) return false;
			if (material != materialType::PEC) return false;
			if (monoDual != MonoDual::Mono) return false;
			if (integralEqu == IntegralEquation::PMCHWT) return false;
			if (movingPathFile.empty()) return false;
			if (!hasMovingRadar) return false;
		}
		if (useFSBMLFMAMoving) {
			if (selectAlgo != selectAlgorithm::MLFMA) return false;
			if (material != materialType::PEC) return false;
			if (monoDual != MonoDual::Mono) return false;
			if (integralEqu == IntegralEquation::PMCHWT) return false;
			if (movingPathFile.empty()) return false;
			if (!hasMovingRadar) return false;
			if (movingAngleMarginDeg < 0.0) return false;
			if (movingAngleMergeTolDeg <= 0.0) return false;
		}
		if (useACARHSMoving) {
			if (selectAlgo != selectAlgorithm::MLFMA) return false;
			if (material != materialType::PEC) return false;
			if (monoDual != MonoDual::Mono) return false;
			if (integralEqu == IntegralEquation::PMCHWT) return false;
			if (movingPathFile.empty()) return false;
			if (!hasMovingRadar) return false;
			if (acaRHSTol <= 0.0) return false;
			if (acaRHSMaxRank < 0) return false;
			if (acaRHSPivotEps <= 0.0) return false;
			if (acaRHSCheckCount < 0) return false;
		}
		return true;
	}

	void validate() const {
		if (!isValid()) {
			throw std::invalid_argument("RCSExportConfig: invalid parameters");
		}
	}

	// ö��ת�ַ����������ļ�����
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

	// �����ļ����ĺ��Ĳ��֣�����·����ǰ׺��ʱ�����
	std::string toCoreString() const {
		std::ostringstream oss;
		oss << algoStr() << "_"
			<< materialStr() << "_"
			<< integralEquStr() << "_"
			<< monoDualStr() << "_"
			<< polWaveStr() << "_"
			<< nPoints << "pts";
		if (useAccurateMLFMAMoving) {
			oss << "_AccurateMLFMAMoving";
		}
		if (useFSBMLFMAMoving) {
			oss << "_FSBMLFMAMoving";
		}
		if (useACARHSMoving) {
			oss << "_ACARHSMoving";
		}
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
	// ----- ͨ���ļ������ɺ��ļ��� -----
	// ��������·���������� + ���Ĵ� + �Զ����׺ + ʱ��� + ��չ��
	std::string buildFilePath(const std::string& suffix, const std::string& extension = "txt") const {
		namespace fs = std::filesystem;
		std::string time = getTimeStamp();
		std::string filename = baseName + "_" + toCoreString() + suffix + "_" + time + "." + extension;
		return (fs::path(parentPath) / filename).string();
	}
	// ֱ�Ӵ��ļ������� ofstream (mode Ĭ�ϲ���Ϊstd::ios::out => �ı�ģʽ|  std::ios::binary => ������ģʽ)
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
