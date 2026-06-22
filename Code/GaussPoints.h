// GaussPoints.h
#pragma once
#ifndef GAUSSPOINT_H
#define GAUSSPOINT_H
class gaussPoints {
public:
	const double* l1;
	const double* l2;
	const double* l3;
	const double* weight;
	gaussPoints(int N_points);
	const int N_points_;
private:
	void selectGauss_Point(const double*& l1, const double*& l2, const double*& l3, const double*& weight, const int N_points);
};
#endif // !GAUSSPOINT_H
