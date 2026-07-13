#pragma once
#ifndef BESSEL_H
#define BESSEL_H

#include <complex>
#include <cmath>
#include <vector>

#include <boost/math/special_functions/bessel.hpp>

namespace Bessel {
	// Spherical Hankel Function of the First Kind
	inline std::complex<double> spherical_hankel_1(int n, double x) {
		double jn = boost::math::sph_bessel(n, x);
		double yn = boost::math::sph_neumann(n, x);
		return std::complex<double>(jn, yn);
	}

	// Spherical Hankel Function of the Second Kind
	inline std::complex<double> spherical_hankel_2(int n, double x) {
		double jn = boost::math::sph_bessel(n, x);
		double yn = boost::math::sph_neumann(n, x);
		return std::complex<double>(jn, -yn);
	}

	// Batch compute Spherical Hankel Function of the Second Kind
	// h_l^(1)(x) = j_l(x) + i*y_l(x), for l = 0, 1, ..., L
	inline void spherical_hankel_1_array(int L, double x,
		std::vector<std::complex<double>>& h1)
	{
		h1.resize(L + 1);

		for (int l = 0; l <= L; ++l) {
			h1[l] = spherical_hankel_1(l, x);
		}
	}

	// h_l^(2)(x) = j_l(x) - i*y_l(x), for l = 0, 1, ..., L
	inline void spherical_hankel_2_array(int L, double x,
		std::vector<std::complex<double>>& h2)
	{
		h2.resize(L + 1);

		for (int l = 0; l <= L; ++l) {
			h2[l] = spherical_hankel_2(l, x);
		}
	}
}// namespace Bessel

namespace ComplexBessel {
	// 复参数球贝塞尔函数 j_n(z) 的初始值
	// j_0(z) = sin(z)/z
	// j_1(z) = sin(z)/z^2 - cos(z)/z
	inline std::complex<double> sph_bessel_j0(const std::complex<double>& z) {
		if (std::abs(z) < 1e-15) return std::complex<double>(1.0, 0.0);
		return std::sin(z) / z;
	}

	inline std::complex<double> sph_bessel_j1(const std::complex<double>& z) {
		if (std::abs(z) < 1e-15) return std::complex<double>(0.0, 0.0);
		return std::sin(z) / (z * z) - std::cos(z) / z;
	}

	// 复参数球诺伊曼函数 y_n(z) 的初始值
	// y_0(z) = -cos(z)/z
	// y_1(z) = -cos(z)/z^2 - sin(z)/z
	inline std::complex<double> sph_neumann_y0(const std::complex<double>& z) {
		return -std::cos(z) / z;
	}

	inline std::complex<double> sph_neumann_y1(const std::complex<double>& z) {
		return -std::cos(z) / (z * z) - std::sin(z) / z;
	}

	// 通过向上递推计算 j_n(z)，对 n <= 某个阈值时稳定
	// 递推关系: f_{n+1}(z) = (2n+1)/z * f_n(z) - f_{n-1}(z)
	// 注意: j_n 的向上递推在 n > |z| 时数值不稳定
	// 建议使用 Miller 算法（向下递推）来计算 j_n

	// Miller 向下递推算法计算 j_0(z) 到 j_N(z)
	inline void sph_bessel_j_array(int N, const std::complex<double>& z,
		std::vector<std::complex<double>>& jn) {
		jn.resize(N + 1);

		if (std::abs(z) < 1e-15) {
			jn[0] = std::complex<double>(1.0, 0.0);
			for (int n = 1; n <= N; ++n) jn[n] = std::complex<double>(0.0, 0.0);
			return;
		}

		// 选择起始点: N_start > N + sqrt(101.5 + N)
		int N_start = N + static_cast<int>(std::sqrt(101.5 + N)) + 10;
		if (N_start < N + 30) N_start = N + 30;

		// 向下递推
		std::complex<double> fn_plus1(0.0, 0.0);
		std::complex<double> fn(1.0, 0.0);  // 任意非零初始值

		std::vector<std::complex<double>> temp(N_start + 1);
		temp[N_start] = fn_plus1;
		temp[N_start - 1] = fn;

		for (int n = N_start - 1; n > 0; --n) {
			std::complex<double> fn_minus1 = (2.0 * n + 1.0) / z * temp[n] - temp[n + 1];
			temp[n - 1] = fn_minus1;
		}

		// 归一化: j_0(z) = sin(z)/z
		std::complex<double> j0_exact = sph_bessel_j0(z);
		std::complex<double> scale = j0_exact / temp[0];

		for (int n = 0; n <= N; ++n) {
			jn[n] = temp[n] * scale;
		}
	}

	// y_n(z) 的向上递推是稳定的
	inline void sph_neumann_y_array(int N, const std::complex<double>& z,
		std::vector<std::complex<double>>& yn) {
		yn.resize(N + 1);
		yn[0] = sph_neumann_y0(z);
		if (N == 0) return;
		yn[1] = sph_neumann_y1(z);
		for (int n = 1; n < N; ++n) {
			yn[n + 1] = (2.0 * n + 1.0) / z * yn[n] - yn[n - 1];
		}
	}

	// 复参数第二类球汉克尔函数 h_n^(2)(z) = j_n(z) - i*y_n(z)
	inline std::complex<double> spherical_hankel_2(int n, const std::complex<double>& z) {
		std::vector<std::complex<double>> jn, yn;
		sph_bessel_j_array(n, z, jn);
		sph_neumann_y_array(n, z, yn);
		return jn[n] - std::complex<double>(0.0, 1.0) * yn[n];
	}

	// 批量计算 h_l^(2)(z) for l = 0, 1, ..., L（更高效）
	inline void spherical_hankel_2_array(int L, const std::complex<double>& z,
		std::vector<std::complex<double>>& h2) {
		std::vector<std::complex<double>> jn, yn;
		sph_bessel_j_array(L, z, jn);
		sph_neumann_y_array(L, z, yn);
		h2.resize(L + 1);
		const std::complex<double> imag_unit(0.0, 1.0);
		for (int l = 0; l <= L; ++l) {
			h2[l] = jn[l] - imag_unit * yn[l];
		}
	}
}

#endif // !BESSEL_H
