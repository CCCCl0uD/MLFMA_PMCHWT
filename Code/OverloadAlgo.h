// OverloadAlgo.h
#pragma once
#ifndef OVERLOADALGO_H
#define OVERLOADALGO_H

#include "PhyConstant.h"
#include "MyStruct.h"
#include "MovingTarget.h"

#include <array>
#include <cmath>
#include <complex>
#include <limits>
#include <set>
#include <utility>
#include <vector>
#include <boost/math/special_functions/legendre.hpp>

// ���غ������ڼ���һ��������صĽ��
/****************************************************************************/
//���
/****************************************************************************/
inline void cross(double* res, double* u, double* v) {
	res[0] = u[1] * v[2] - u[2] * v[1];
	res[1] = u[2] * v[0] - u[0] * v[2];
	res[2] = u[0] * v[1] - u[1] * v[0];
}

inline std::array<double, 3> cross(
	const std::array<double, 3>& u,
	const std::array<double, 3>& v)
{
	return {
		u[1] * v[2] - u[2] * v[1],
		u[2] * v[0] - u[0] * v[2],
		u[0] * v[1] - u[1] * v[0]
	};
}

inline void cross(double* res, const double* u, const double* v) {
	res[0] = u[1] * v[2] - u[2] * v[1];
	res[1] = u[2] * v[0] - u[0] * v[2];
	res[2] = u[0] * v[1] - u[1] * v[0];
}

inline void cross(std::complex<double>* res, double* u, std::complex<double>* v)
{
	res[0] = u[1] * v[2] - u[2] * v[1];
	res[1] = u[2] * v[0] - u[0] * v[2];
	res[2] = u[0] * v[1] - u[1] * v[0];
}

inline void cross(std::complex<double>* res, std::complex<double>* u, std::complex<double>* v)
{
	res[0] = u[1] * v[2] - u[2] * v[1];
	res[1] = u[2] * v[0] - u[0] * v[2];
	res[2] = u[0] * v[1] - u[1] * v[0];
}

inline void cross(Complex3D& res, double* u, Complex3D& v) {
	res.x = u[1] * v.z - u[2] * v.y;
	res.y = u[2] * v.x - u[0] * v.z;
	res.z = u[0] * v.y - u[1] * v.x;
}

inline void cross(Complex3D& res, Complex3D& u, Complex3D& v) {
	res.x = u.y * v.z - u.z * v.y;
	res.y = u.z * v.x - u.x * v.z;
	res.z = u.x * v.y - u.y * v.x;
}

inline void cross(Complex3D& res, const double* u, const Complex3D& v) {
	res.x = u[1] * v.z - u[2] * v.y;
	res.y = u[2] * v.x - u[0] * v.z;
	res.z = u[0] * v.y - u[1] * v.x;
}

/****************************************************************************/
//���
/****************************************************************************/
template<typename T, typename C>
inline C dot(T* u, C* v) {
	return u[0] * v[0] + u[1] * v[1] + u[2] * v[2];
}

template<typename T, typename C>
inline T dot(T* u, C& v) {
	return u[0] * v[0] + u[1] * v[1] + u[2] * v[2];
}

inline std::complex<double> dot(const Complex3D& u, const Complex3D& v) {
	return u.x * v.x + u.y * v.y + u.z * v.z;
}

template<typename C>
inline C dot_conj(Complex3D& u, Complex3D& v) {
	return u.x * std::conj(v.x) + u.y * std::conj(v.y) + u.z * std::conj(v.z);
}
/****************************************************************************/
//����
/****************************************************************************/
template<typename T>
inline void div(T* res, T* u, double& CONST) {
	res[0] = u[0] / CONST;
	res[1] = u[1] / CONST;
	res[2] = u[2] / CONST;
}

template<typename T>
inline void div(T& res, T& u, double& CONST) {
	res[0] = u[0] / CONST;
	res[1] = u[1] / CONST;
	res[2] = u[2] / CONST;
}
/****************************************************************************/
//�˷�
/****************************************************************************/
template<typename T, typename C>
inline void mul(T* res, C* u, T& Const) {
	res[0] = u[0] * Const;
	res[1] = u[1] * Const;
	res[2] = u[2] * Const;
}

template<typename T>
inline void mul(T* res, T* u, T& Const) {
	res[0] = u[0] * Const;
	res[1] = u[1] * Const;
	res[2] = u[2] * Const;
}

template<typename T>
inline void mul(T& res, T& u, double& Const) {
	res[0] = u[0] * Const;
	res[1] = u[1] * Const;
	res[2] = u[2] * Const;
}

template<typename T, typename C>
inline void mul(C* res, C* u, T& Const) {
	res[0] = u[0] * Const;
	res[1] = u[1] * Const;
	res[2] = u[2] * Const;
}

template<typename T, typename C>
inline void mul(C* res, T* u, T& Const) {
	res[0] = u[0] * Const;
	res[1] = u[1] * Const;
	res[2] = u[2] * Const;
}

/****************************************************************************/
//�ӷ�
/****************************************************************************/
template<typename T, typename C>
inline void add(C* res, T* u, C* v)
{
	res[0] = u[0] + v[0];
	res[1] = u[1] + v[1];
	res[2] = u[2] + v[2];
}

template<typename T, typename C>
inline void add(C& res, T& u, C& v) {
	res[0] = u[0] + v[0];
	res[1] = u[1] + v[1];
	res[2] = u[2] + v[2];
}

template<typename C>
inline void add_cp(C* res, C* u, C* v) {
	res[0] = (u[0] + v[0]) / 2;
	res[1] = (u[1] + v[1]) / 2;
	res[2] = (u[2] + v[2]) / 2;
}
/****************************************************************************/
//����
/****************************************************************************/
template<typename T>
inline void sub(T* res, T* u, T* v)
{
	res[0] = u[0] - v[0];
	res[1] = u[1] - v[1];
	res[2] = u[2] - v[2];
}

template<typename T>
inline void sub(T& res, T& u, T& v) {
	res[0] = u[0] - v[0];
	res[1] = u[1] - v[1];
	res[2] = u[2] - v[2];
}
/****************************************************************************/
//ģ��
/****************************************************************************/
template<typename T>
inline T norm2(T* vec) {
	return dot<T>(vec, vec);
}

template<typename T>
inline T norm(T* vec) {
	return sqrt(norm2<T>(vec));
}

template<typename T>
inline T norm2(T* u, T* v)
{
	const T tmp0 = u[0] - v[0];
	const T tmp1 = u[1] - v[1];
	const T tmp2 = u[2] - v[2];
	return tmp0 * tmp0 + tmp1 * tmp1 + tmp2 * tmp2;
}

template<typename T>
inline  T norm(T* u, T* v)
{
	return sqrt(norm2(u, v));
}

namespace my_memory {
	template<typename T>
	inline size_t memory4D(const std::vector<std::vector<std::vector<std::vector<T>>>>& V)
	{
		size_t mem = sizeof(V);
		for (const auto& l1 : V)
			for (const auto& l2 : l1)
				for (const auto& l3 : l2)
					mem += sizeof(T) * l3.capacity();
		return mem;
	}

	template<typename T>
	inline size_t memory3D(const std::vector<std::vector<std::vector<T>>>& V) {
		size_t mem = sizeof(V);
		for (const auto& lvl : V)
			for (const auto& row : lvl)
				mem += sizeof(T) * row.capacity();
		return mem;
	}

	template<typename T>
	inline size_t memory2D(const std::vector<std::vector<T>>& V) {
		size_t mem = sizeof(V);
		for (const auto& row : V) {
			mem += sizeof(T) * row.capacity();
		}
		return mem;
	}

	template<typename T>
	inline size_t memory1D(const std::vector<T>& V)
	{
		return sizeof(V) + sizeof(T) * V.capacity();
	}
}// namespace my_memory

namespace my_math {
	// ţ��-����ѷ������Ѱ�����õ¶���ʽ�ĸ�
	inline double newton_raphson(int n, double x0, double tolerance, int max_iter) {
		double x = x0;
		for (int i = 0; i < max_iter; ++i) {
			double p = boost::math::legendre_p(n, x);
			double dp = boost::math::legendre_p_prime(n, x);
			double dx = p / dp;
			x -= dx;
			if (std::abs(dx) < tolerance) {
				break;
			}
		}
		return x;
	}

	// ����n���˹-���õ»��ֵĽڵ��Ȩ��
	inline void gauss_legendre(int n, std::vector<double>& XGL_temp, std::vector<double>& WGL_temp) {
		XGL_temp.resize(n);
		WGL_temp.resize(n);
		const double tolerance = 1e-10;
		const int max_iter = 100;
		for (int i = 0; i < (n + 1) / 2; ++i) {
			// ��ʼ�²�ֵ
			double x0 = std::cos(Pi * (i + 0.75) / (n + 0.5));
			double root = newton_raphson(n, x0, tolerance, max_iter);
			XGL_temp[i] = -root;
			XGL_temp[n - i - 1] = root;

			double dp = boost::math::legendre_p_prime(n, root);
			double weight = 2.0 / ((1 - root * root) * dp * dp);
			WGL_temp[i] = weight;
			WGL_temp[n - i - 1] = weight;
		}
	}

	// ���õ¶���ʽ�������ڻ���
	inline double legendrePolynomialParameters(double r[3], double kp[3]) {
		return r[0] * kp[0] + r[1] * kp[1] + r[2] * kp[2];
	}

	inline double lagrangeWeight(double x0, double x1, double x2, double x3) {
		return ((x0 - x2) * (x0 - x3)) / ((x1 - x2) * (x1 - x3));
	}

	// ���� (n �� h) �� r
	template<typename Vec1, typename Vec2, typename Vec3>
	inline double cross_dot(const Vec1& n, const Vec2& h, const Vec3& r) {
		// ����������
		double cx = n[1] * h[2] - n[2] * h[1];
		double cy = n[2] * h[0] - n[0] * h[2];
		double cz = n[0] * h[1] - n[1] * h[0];
		// ���
		return cx * r[0] + cy * r[1] + cz * r[2];
	}

	static const std::complex<double> jl_table[4] = {
	{1.0, 0.0}, {0.0, -1.0}, {-1.0, 0.0}, {0.0, 1.0}
	};

	inline void normalize(double r[3], double length) {
		r[0] /= length;
		r[1] /= length;
		r[2] /= length;
	}

	inline int calculateNumPoints(const double sca_th_s, const double sca_ph_s, const double sca_th_f, const double sca_ph_f, const double step) {
		// ����th����ĵ�����ph����ĵ���
		int num_th_points, num_ph_points;

		if (fabs(sca_th_f - sca_th_s) > 1e-10) {
			// th����仯������th����ĵ���
			num_th_points = static_cast<int>(fabs(sca_th_f - sca_th_s) / step) + 1;
			// ph���򲻱䣬ֻ��1����
			num_ph_points = 1;
		}
		else if (fabs(sca_ph_f - sca_ph_s) > 1e-10) {
			// ph����仯������ph����ĵ���
			num_ph_points = static_cast<int>(fabs(sca_ph_f - sca_ph_s) / step) + 1;
			// th���򲻱䣬ֻ��1����
			num_th_points = 1;
		}
		else {
			// �������򶼲��䣬ֻ��1����
			num_th_points = 1;
			num_ph_points = 1;
		}
		// �ܵ��� = th������� �� ph�������
		return num_th_points * num_ph_points;
	}

	inline void computeTensorDotVector(std::array<std::complex<double>, 3>& F, const std::array<double, 3>& r, const std::array<double, 3>& rho) {
		const double r0r0 = r[0] * r[0];
		const double r0r1 = r[0] * r[1];
		const double r0r2 = r[0] * r[2];
		const double r1r0 = r[1] * r[0];
		const double r1r1 = r[1] * r[1];
		const double r1r2 = r[1] * r[2];
		const double r2r0 = r[2] * r[0];
		const double r2r1 = r[2] * r[1];
		const double r2r2 = r[2] * r[2];

		F[0] = (1 - r0r0) * rho[0] - r0r1 * rho[1] - r0r2 * rho[2];
		F[1] = -r1r0 * rho[0] + (1 - r1r1) * rho[1] - r1r2 * rho[2];
		F[2] = -r2r0 * rho[0] - r2r1 * rho[1] + (1 - r2r2) * rho[2];
	}
}// namespace my_math

namespace my_fsb {
	inline double normalizePhiDeg(double phiDeg) {
		double v = std::fmod(phiDeg, 360.0);
		if (v < 0.0) v += 360.0;
		return v;
	}

	inline double wrapAngleRad(double x) {
		while (x > Pi) x -= 2.0 * Pi;
		while (x < -Pi) x += 2.0 * Pi;
		return x;
	}

	inline double dirichletKernel(double delta, int order) {
		const int n = std::max(1, order);
		const double half = 0.5 * delta;
		const double denom = std::sin(half);
		if (std::abs(denom) < 1.0e-12) {
			return 1.0;
		}
		return std::sin((2.0 * n + 1.0) * half) /
			((2.0 * n + 1.0) * denom);
	}

	template<typename SolverType>
	inline double estimateBoundingSphereRadius(const SolverType& solver) {
		std::array<double, 3> mn = {
			std::numeric_limits<double>::max(),
			std::numeric_limits<double>::max(),
			std::numeric_limits<double>::max()
		};
		std::array<double, 3> mx = {
			std::numeric_limits<double>::lowest(),
			std::numeric_limits<double>::lowest(),
			std::numeric_limits<double>::lowest()
		};

		auto visitPoint = [&](const Point* p) {
			if (!p) return;
			mn[0] = std::min(mn[0], p->x);
			mn[1] = std::min(mn[1], p->y);
			mn[2] = std::min(mn[2], p->z);
			mx[0] = std::max(mx[0], p->x);
			mx[1] = std::max(mx[1], p->y);
			mx[2] = std::max(mx[2], p->z);
			};

		for (const auto& rwg : solver.rwgs) {
			if (rwg.positiveFace) {
				visitPoint(rwg.positiveFace->vertex1);
				visitPoint(rwg.positiveFace->vertex2);
				visitPoint(rwg.positiveFace->vertex3);
			}
			if (rwg.negativeFace) {
				visitPoint(rwg.negativeFace->vertex1);
				visitPoint(rwg.negativeFace->vertex2);
				visitPoint(rwg.negativeFace->vertex3);
			}
		}

		if (mx[0] < mn[0] || mx[1] < mn[1] || mx[2] < mn[2]) {
			return 0.0;
		}

		const std::array<double, 3> center = {
			0.5 * (mn[0] + mx[0]),
			0.5 * (mn[1] + mx[1]),
			0.5 * (mn[2] + mx[2])
		};

		double radius = 0.0;
		auto updateRadius = [&](const Point* p) {
			if (!p) return;
			const double dx = p->x - center[0];
			const double dy = p->y - center[1];
			const double dz = p->z - center[2];
			radius = std::max(radius, std::sqrt(dx * dx + dy * dy + dz * dz));
			};

		for (const auto& rwg : solver.rwgs) {
			if (rwg.positiveFace) {
				updateRadius(rwg.positiveFace->vertex1);
				updateRadius(rwg.positiveFace->vertex2);
				updateRadius(rwg.positiveFace->vertex3);
			}
			if (rwg.negativeFace) {
				updateRadius(rwg.negativeFace->vertex1);
				updateRadius(rwg.negativeFace->vertex2);
				updateRadius(rwg.negativeFace->vertex3);
			}
		}
		return radius;
	}

	inline double deltaPhiDeg(double a, double b) {
		const double da = normalizePhiDeg(a);
		const double db = normalizePhiDeg(b);
		const double d = std::abs(da - db);
		return std::min(d, 360.0 - d);
	}

	inline double movingThetaRingDeg(int ringIndex, int thetaCount) {
		return 180.0 * static_cast<double>(ringIndex + 1) /
			static_cast<double>(thetaCount + 1);
	}

	inline double movingPhiSampleDeg(int phiIndex, int phiOrder) {
		const int phiCount = 2 * phiOrder + 1;
		return normalizePhiDeg(
			360.0 * static_cast<double>(phiIndex) /
			static_cast<double>(phiCount));
	}

	inline std::vector<int> selectMovingThetaRingIndices(
		double thetaQueryDeg,
		int thetaCount,
		double marginDeg)
	{
		std::vector<int> indices;
		if (thetaCount <= 0) {
			return indices;
		}

		const double spacingDeg = 180.0 / static_cast<double>(thetaCount + 1);
		const double radiusDeg = std::max(0.0, marginDeg) + 0.5 * spacingDeg;
		int nearest = 0;
		double nearestDiff = std::numeric_limits<double>::max();
		for (int i = 0; i < thetaCount; ++i) {
			const double thetaDeg = movingThetaRingDeg(i, thetaCount);
			const double diff = std::abs(thetaDeg - thetaQueryDeg);
			if (diff <= radiusDeg) {
				indices.push_back(i);
			}
			if (diff < nearestDiff) {
				nearestDiff = diff;
				nearest = i;
			}
		}
		if (indices.empty()) {
			indices.push_back(nearest);
		}
		return indices;
	}

	inline std::vector<int> selectMovingPhiIndices(
		double phiQueryDeg,
		int phiOrder,
		double marginDeg)
	{
		std::vector<int> indices;
		const int safeOrder = std::max(1, phiOrder);
		const int phiCount = 2 * safeOrder + 1;
		const double spacingDeg = 360.0 / static_cast<double>(phiCount);
		const double radiusDeg = std::max(0.0, marginDeg) + 0.5 * spacingDeg;
		int nearest = -safeOrder;
		double nearestDiff = std::numeric_limits<double>::max();

		for (int n = -safeOrder; n <= safeOrder; ++n) {
			const double phiDeg = movingPhiSampleDeg(n, safeOrder);
			const double diff = deltaPhiDeg(phiDeg, phiQueryDeg);
			if (diff <= radiusDeg) {
				indices.push_back(n);
			}
			if (diff < nearestDiff) {
				nearestDiff = diff;
				nearest = n;
			}
		}
		if (indices.empty()) {
			indices.push_back(nearest);
		}
		return indices;
	}

	inline std::vector<std::vector<FSBMonoSample>> buildMovingFSBSampleLayout(
		const std::vector<MovingAngleSample>& queryAngles,
		double W,
		double marginDeg)
	{
		const int thetaCount = std::max(2, static_cast<int>(std::ceil(2.0 * std::max(1.0, W))));
		std::vector<std::set<int>> selectedPhi(thetaCount);

		for (const auto& query : queryAngles) {
			const std::vector<int> thetaIndices = selectMovingThetaRingIndices(
				query.thetaBodyDeg, thetaCount, marginDeg);
			for (int thetaIndex : thetaIndices) {
				const double theta = movingThetaRingDeg(thetaIndex, thetaCount) * Pi / 180.0;
				const int phiOrder = std::max(
					1, static_cast<int>(std::ceil(std::max(1.0, W) * std::sin(theta))));
				const std::vector<int> phiIndices = selectMovingPhiIndices(
					query.phiBodyDeg, phiOrder, marginDeg);
				selectedPhi[thetaIndex].insert(phiIndices.begin(), phiIndices.end());
			}
		}

		std::vector<std::vector<FSBMonoSample>> rings;
		rings.reserve(thetaCount);
		for (int thetaIndex = 0; thetaIndex < thetaCount; ++thetaIndex) {
			const double thetaDeg = movingThetaRingDeg(thetaIndex, thetaCount);
			const double theta = thetaDeg * Pi / 180.0;
			const int phiOrder = std::max(
				1, static_cast<int>(std::ceil(std::max(1.0, W) * std::sin(theta))));
			std::vector<FSBMonoSample> ring;
			ring.reserve(selectedPhi[thetaIndex].size());
			for (int phiIndex : selectedPhi[thetaIndex]) {
				FSBMonoSample sample;
				sample.thetaDeg = thetaDeg;
				sample.phiDeg = movingPhiSampleDeg(phiIndex, phiOrder);
				sample.phiOrder = phiOrder;
				ring.push_back(std::move(sample));
			}
			rings.push_back(std::move(ring));
		}
		return rings;
	}
	inline std::complex<double> interpolateFSBComponent(
		const std::vector<std::vector<FSBMonoSample>>& rings,
		double thetaDeg,
		double phiDeg,
		int component)
	{
		const double theta = thetaDeg * Pi / 180.0;
		const double phi = normalizePhiDeg(phiDeg) * Pi / 180.0;
		const int thetaOrder = std::max(1, static_cast<int>(rings.size()));

		std::complex<double> result(0.0, 0.0);
		double thetaWeightSum = 0.0;

		for (const auto& ring : rings) {
			if (ring.empty()) continue;
			const double thetaM = ring.front().thetaDeg * Pi / 180.0;
			const double wt = dirichletKernel(theta - thetaM, thetaOrder);
			std::complex<double> ringValue(0.0, 0.0);
			double phiWeightSum = 0.0;

			for (const auto& sample : ring) {
				const double phiN = normalizePhiDeg(sample.phiDeg) * Pi / 180.0;
				const double wp = dirichletKernel(
					wrapAngleRad(phi - phiN), sample.phiOrder);
				ringValue += wp * sample.field[component];
				phiWeightSum += wp;
			}

			if (std::abs(phiWeightSum) > 1.0e-12) {
				ringValue /= phiWeightSum;
			}
			result += wt * ringValue;
			thetaWeightSum += wt;
		}

		if (std::abs(thetaWeightSum) > 1.0e-12) {
			result /= thetaWeightSum;
		}
		return result;
	}

	inline std::vector<std::complex<double>> interpolateFSBField(
		const std::vector<std::vector<FSBMonoSample>>& rings,
		double thetaDeg,
		double phiDeg)
	{
		std::vector<std::complex<double>> field(3, std::complex<double>(0.0, 0.0));
		for (int c = 0; c < 3; ++c) {
			field[c] = interpolateFSBComponent(rings, thetaDeg, phiDeg, c);
		}
		return field;
	}
}

#endif // OVERLOADALGO_H