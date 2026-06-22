#include "EMSource.h"

EMSource::EMSource() : EMSource(3.0e8, 0.0) {
}

EMSource::EMSource(double freq, double pol) {
	data_fre = freq;
	data_pol = pol;
	data_waveLength = c0 / freq;
	data_k = 2.0 * Pi / data_waveLength;
	data_k1.real(data_k);   data_k1.imag(0.0);
	data_eta1.real(eta0_);   data_eta1.imag(0.0);
	data_pol = change(pol);
}

EMSource::~EMSource() {
}

void EMSource::initDielectric(
	const std::complex<double>& epsR,
	const std::complex<double>& muR)
{
	// ÓÐºÄ½éÖÊ £¨¦År = ¦Å' - j¦Å"£©
	data_isDielectric = true;
	data_epsilonR = epsR;
	data_muR = muR;
	data_k1 = std::complex<double>(data_k, 0.0);
	data_k2 = std::complex<double>(data_k, 0.0) * std::sqrt(epsR * muR);
	data_etai = std::sqrt(muR / epsR);
	data_eta1 = std::complex<double>(eta0_, 0.0);
	data_eta2 = std::complex<double>(eta0_, 0.0) * std::sqrt(muR / epsR);
}

void EMSource::initPW(double pol, double th, double ph) {
	double u[3], v[3];
	data_pol = change(pol);
	data_th = change(th);
	data_ph = change(ph);
	data_vW[0] = -std::sin(data_th) * std::cos(data_ph);//»¡¶È
	data_vW[1] = -std::sin(data_th) * std::sin(data_ph);
	data_vW[2] = -std::cos(data_th);
	returnZero(data_vW[0]);
	returnZero(data_vW[1]);
	returnZero(data_vW[2]);

	data_vTh[0] = -std::cos(data_th) * std::cos(data_ph);
	data_vTh[1] = -std::cos(data_th) * std::sin(data_ph);
	data_vTh[2] = std::sin(data_th);
	returnZero(data_vTh[0]);
	returnZero(data_vTh[1]);
	returnZero(data_vTh[2]);

	data_vPh[0] = -std::sin(data_ph);
	data_vPh[1] = std::cos(data_ph);
	data_vPh[2] = 0.0;
	returnZero(data_vPh[0]);
	returnZero(data_vPh[1]);
	returnZero(data_vPh[2]);
	if (compare(data_pol, change(90.0))) {
		data_vE[0] = data_vPh[0];
		data_vE[1] = data_vPh[1];
		data_vE[2] = data_vPh[2];
		crossPW(data_vW, data_vE, data_vH);
		normalize(data_vH);
	}
	if (compare(data_pol, 0.0)) {
		data_vE[0] = data_vTh[0];
		data_vE[1] = data_vTh[1];
		data_vE[2] = data_vTh[2];
		crossPW(data_vW, data_vE, data_vH);
		normalize(data_vH);
	}
}