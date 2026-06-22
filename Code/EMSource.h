#pragma once
#ifndef EMSOURCE_H
#define EMSOURCE_H
#include "PhyConstant.h"

class EMSource
{
public:
	EMSource();
	EMSource(double freq, double pol);
	~EMSource();

	double wavelength() const { return data_waveLength; }
	double freq() const { return data_fre; }

	void initDielectric(const std::complex<double>& epsR, const std::complex<double>& muR);
	bool isDielectric() const { return data_isDielectric; }
	std::complex<double> k1() const { return data_k1; }
	std::complex<double> k2() const { return data_k2; }
	std::complex<double> eta1() const { return data_eta1; }
	std::complex<double> eta2() const { return data_eta2; }
	std::complex<double> etai() const { return data_etai; }
	std::complex<double> epsilonR() const { return data_epsilonR; }
	std::complex<double> muR() const { return data_muR; }

	double k1_abs() const { return std::abs(data_k1); }
	double k1_real() const { return data_k1.real(); }
	double k2_abs() const { return std::abs(data_k2); }
	double k2_real() const { return data_k2.real(); }

	double& pol() { return data_pol; }
	double& th() { return data_th; }
	double& ph() { return data_ph; }
	double& vTh(int coordinate) { return data_vTh[coordinate]; }
	double& vPh(int coordinate) { return data_vPh[coordinate]; }
	double& vW(int coordinate) { return data_vW[coordinate]; }
	double& vE(int coordinate) { return data_vE[coordinate]; }
	double& vH(int coordinate) { return data_vH[coordinate]; }

	void initPW(double pol, double th, double ph);

	double change(double a) {
		double b;
		b = a * Pi / 180.0;
		if (std::abs(b) < 1e-6) {
			b = 0.0;
		}
		return b;
	}

	bool compare(double a, double b) {
		if (std::abs(a - b) < 1e-6) {
			return true;
		}
		return false;
	}

	void returnZero(double& a) {
		if (std::abs(a) < 1e-6) {
			a = 0.0;
		}
	}

	inline void crossPW(const double a[3], const double b[3], double c[3]) {
		c[0] = a[1] * b[2] - a[2] * b[1];
		c[1] = a[2] * b[0] - a[0] * b[2];
		c[2] = a[0] * b[1] - a[1] * b[0];
	}

	inline void normalize(double a[3]) {
		double norm = std::sqrt(a[0] * a[0] + a[1] * a[1] + a[2] * a[2]);
		if (norm > 1e-15) {
			a[0] /= norm;
			a[1] /= norm;
			a[2] /= norm;
		}
	}

protected:
	const double epsilon0 = 8.8541878128e-12;// 真空介电常数 (F/m)
	const double mu0 = 1.2566370614359173e-6;// 真空磁导率 (H/m)
	const double eta0_ = 119.9169832 * Pi;	 // 真空中波阻抗 (Ω)

	double data_fre, data_pol, data_waveLength, data_k;
	double data_th, data_ph;
	double data_vTh[3], data_vPh[3];
	double data_vW[3], data_vE[3], data_vH[3];

	bool data_isDielectric = false;
	std::complex<double> data_k1 = std::complex<double>(0.0, 0.0);
	std::complex<double> data_k2 = std::complex<double>(0.0, 0.0);
	std::complex<double> data_eta1 = std::complex<double>(eta0_, 0.0);
	std::complex<double> data_eta2 = std::complex<double>(eta0_, 0.0);
	std::complex<double> data_etai = std::complex<double>(1.0, 0.0);
	std::complex<double> data_epsilonR = std::complex<double>(1.0, 0.0);
	std::complex<double> data_muR = std::complex<double>(1.0, 0.0);
};

#endif // EMSOURCE_H