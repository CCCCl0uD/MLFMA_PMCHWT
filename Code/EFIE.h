// EFIE.h
#pragma once
#ifndef EFIE_H
#define EFIE_H

#include "MyStruct.h"
#include "GaussPoints.h"

namespace EFIE {
	void computeEFIE_Zij(const RWGBase* rwgField, const RWGBase* rwgSource, const std::complex<double> k_, const std::complex<double> eta_,
		const gaussPoints& gp, std::complex<double>& Zij_near_E);

	void singularityEFIE(double& areaTriF, double& areaTriS,
		double* vertexF1, double* vertexF2, double* vertexF3, double* nonPublicEdge_vertexF, const double recAf,
		double* vertexS1, double* vertexS2, double* vertexS3, double* nonPublicEdge_vertexS, const double recAs, double* normalVectorS,
		const gaussPoints& gp, const std::complex<double> k_, const std::complex<double> eta_, std::complex<double>& alok);
}
#endif // !EFIE_H
