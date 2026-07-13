// PhyConstant.h
#pragma once
#ifndef PHYCONSTANT_H
#define PHYCONSTANT_H
#include <complex>

const double Pi = 3.14159265358979323846;// 圆周率

const double c0 = 299792458.0;			 // 真空中的光速 (m/s)

const std::complex<double> J(0.0, 1.0);	 // 虚数单位

const double eps0 = 2.2204e-16;			 // 浮点数截断数
const double leafCoeff_ = 0.3;           // 最小盒子边长与波长的关系系数
const double alpha = 0.2;				 // EFIE和CFIE组合时的权重因子

#endif // PHYCONSTANT_H
