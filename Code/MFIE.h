// MFIE.h
#pragma once
#ifndef MFIE_H
#define MFIE_H

#include "MyStruct.h"
#include "GaussPoints.h"

namespace MFIE {
	void computeMFIE_Zij(const RWGBase* rwgField, const RWGBase* rwgSource, const std::complex<double> k_,
		const gaussPoints& gp, std::complex<double>& Zij_near_M);

	void singularityMFIE(double& areaTriF, double& areaTriS,
		double* vertexF1, double* vertexF2, double* vertexF3,
		double* nonPublicEdge_vertexF, double* normalVectorF, const double recAf,
		double* vertexS1, double* vertexS2, double* vertexS3,
		double* nonPublicEdge_vertexS, double* normalVectorS, const double recAs,
		const gaussPoints& gp, const std::complex<double> k_, std::complex<double>& alok);
}
#endif // !MFIE_H
