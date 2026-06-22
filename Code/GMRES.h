// MGMRES.h
#pragma once
#ifndef GMRES_H
#define GMRES_H

#include <algorithm>
#include <cmath>
#include <complex>
#include <vector>
#include <iostream>

namespace MGMRES {
	// 辅助函数：共轭点积
	inline std::complex<double> dot_conj(int n, const std::complex<double>* a, const std::complex<double>* b) {
		std::complex<double> sum = 0.0;
		for (int i = 0; i < n; ++i)
			sum += std::conj(a[i]) * b[i];
		return sum;
	}

	// 生成 Givens 旋转
	inline void generate_givens(std::complex<double> a, std::complex<double> b,
		std::complex<double>& c, std::complex<double>& s) {
		const double eps = 1e-16;
		double a_norm = std::abs(a);
		double b_norm = std::abs(b);
		if (b_norm < eps) {
			c = 1.0;
			s = 0.0;
		}
		else if (a_norm < eps) {
			c = 0.0;
			s = 1.0;
		}
		else {
			double r = std::sqrt(std::real(a * std::conj(a) + b * std::conj(b)));
			c = a / r;
			s = b / r;
		}
	}

	// 应用 Givens 旋转
	inline void apply_givens(std::complex<double>& x, std::complex<double>& y,
		std::complex<double> c, std::complex<double> s) {
		std::complex<double> temp = std::conj(c) * x + std::conj(s) * y;
		y = -s * x + c * y;
		x = temp;
	}

	/**
	* 广义最小残差法 (GMRES) 迭代求解器
	*
	* @tparam AxFunc 矩阵向量乘法函数类型，签名：void (int n, const std::complex<double>* x, std::complex<double>* w)
	* @param n       矩阵维度
	* @param x       解向量（输入输出）
	* @param rhs     右端项
	* @param itr_max 最大迭代次数
	* @param mr      Krylov 子空间重启大小
	* @param tol_abs 绝对容差
	* @param tol_rel 相对容差
	* @param ax      矩阵向量乘法函数对象
	* @param verbose 是否输出迭代信息
	* @return        true 表示收敛，false 表示达到最大迭代次数
	*/
	template<typename AxFunc>
	bool mgmres(int n, std::complex<double>* x, const std::complex<double>* rhs, int itr_max, int mr, double tol_abs, double tol_rel, AxFunc ax, bool verbose = true) {
		using Complex = std::complex<double>;
		const double eps = 1e-16;

		// 工作向量
		std::vector<Complex> r(n), w(n), z(n);
		std::vector<Complex> v(n * (mr + 1));
		std::vector<Complex> h((mr + 1) * mr, Complex(0.0, 0.0));
		std::vector<Complex> cs(mr), sn(mr), s(mr + 1);
		std::vector<Complex> y(mr);

		int total_iter = 0;
		double error = 0.0;

		// 初始残差 r = rhs - A*x
		ax(n, x, w.data());                 // w = A*x
		for (int i = 0; i < n; ++i)
			r[i] = rhs[i] - w[i];

		double bnorm = std::sqrt(std::real(dot_conj(n, rhs, rhs)));
		if (bnorm < eps) bnorm = 1.0;

		double rnorm = std::sqrt(std::real(dot_conj(n, r.data(), r.data())));
		error = rnorm / bnorm;
		if (error < tol_rel || error < tol_abs) return true;

		// 外层重启循环
		for (int restart_cycle = 1; restart_cycle <= itr_max; ++restart_cycle) {
			// Arnoldi 初始化
			for (int i = 0; i < n; ++i) v[i] = r[i] / rnorm;
			std::fill(s.begin(), s.end(), Complex(0.0, 0.0));
			s[0] = rnorm;

			// 重置 Hessenberg 矩阵和旋转系数
			std::fill(h.begin(), h.end(), Complex(0.0, 0.0));
			std::fill(cs.begin(), cs.end(), Complex(0.0, 0.0));
			std::fill(sn.begin(), sn.end(), Complex(0.0, 0.0));

			int k_used = 0;
			bool converged_inner = false;

			// 内层 Krylov 子空间构造
			for (int k = 0; k < mr; ++k) {
				if (total_iter >= itr_max) break;
				total_iter++;
				k_used = k + 1;

				// z = v_k
				for (int i = 0; i < n; ++i) z[i] = v[i + k * n];

				// 矩阵向量乘法 w = A * z
				ax(n, z.data(), w.data());

				// 修正 Gram-Schmidt 带重新正交化
				for (int j = 0; j <= k; ++j) {
					h[j + k * (mr + 1)] = dot_conj(n, &v[j * n], w.data());
					for (int i = 0; i < n; ++i)
						w[i] -= h[j + k * (mr + 1)] * v[i + j * n];
				}
				// 重新正交化
				for (int j = 0; j <= k; ++j) {
					Complex temp = dot_conj(n, &v[j * n], w.data());
					h[j + k * (mr + 1)] += temp;
					for (int i = 0; i < n; ++i)
						w[i] -= temp * v[i + j * n];
				}

				double h_next = std::sqrt(std::real(dot_conj(n, w.data(), w.data())));
				h[k + 1 + k * (mr + 1)] = h_next;

				if (h_next < eps * rnorm) {
					if (verbose) std::cout << "MGMRES warning: Krylov subspace linearly dependent\n";
					break;
				}

				for (int i = 0; i < n; ++i)
					v[i + (k + 1) * n] = w[i] / h_next;

				// 应用已有的 Givens 旋转
				for (int j = 0; j < k; ++j) {
					apply_givens(h[j + k * (mr + 1)], h[j + 1 + k * (mr + 1)], cs[j], sn[j]);
				}
				// 生成新的 Givens 旋转并应用
				generate_givens(h[k + k * (mr + 1)], h[k + 1 + k * (mr + 1)], cs[k], sn[k]);
				apply_givens(h[k + k * (mr + 1)], h[k + 1 + k * (mr + 1)], cs[k], sn[k]);
				h[k + 1 + k * (mr + 1)] = Complex(0.0, 0.0);
				apply_givens(s[k], s[k + 1], cs[k], sn[k]);

				error = std::abs(s[k + 1]) / bnorm;
				if (verbose) {
					std::cout << "Cycle " << restart_cycle << "." << k + 1 << " Res: " << error << "\n";
				}
				if (error < tol_rel || error < tol_abs) {
					converged_inner = true;
					break;
				}
			}

			// 回代求解最小二乘问题
			int m = std::min(k_used, mr);
			std::vector<Complex> h_temp(m * m, Complex(0.0, 0.0));
			std::vector<Complex> s_temp(m, Complex(0.0, 0.0));
			for (int i = 0; i < m; ++i) {
				for (int j = i; j < m; ++j)
					h_temp[i + j * m] = h[i + j * (mr + 1)];
				s_temp[i] = s[i];
			}
			// 上三角回代
			for (int i = m - 1; i >= 0; --i) {
				y[i] = s_temp[i];
				for (int j = i + 1; j < m; ++j)
					y[i] -= h_temp[i + j * m] * y[j];
				if (std::abs(h_temp[i + i * m]) > eps)
					y[i] /= h_temp[i + i * m];
				else
					y[i] = Complex(0.0, 0.0);
			}

			// 更新解 x = x + V * y
			for (int j = 0; j < m; ++j) {
				for (int i = 0; i < n; ++i)
					x[i] += y[j] * v[i + j * n];
			}

			// 精确残差用于重启
			ax(n, x, w.data());
			for (int i = 0; i < n; ++i) r[i] = rhs[i] - w[i];
			rnorm = std::sqrt(std::real(dot_conj(n, r.data(), r.data())));
			error = rnorm / bnorm;

			if (verbose && !converged_inner) {
				std::cout << "Restart cycle " << restart_cycle << " exact residual: " << error << "\n";
			}

			if (error < tol_rel || error < tol_abs) {
				if (verbose) std::cout << "GMRES converged! Final residual: " << error << "\n";
				return true;
			}
			if (total_iter >= itr_max) {
				if (verbose) std::cout << "GMRES reached max iterations, final residual: " << error << "\n";
				return false;
			}
			if (error > 100.0 * tol_rel && restart_cycle > 3) {
				if (verbose) std::cout << "Warning: residual growing rapidly, algorithm may be unstable\n";
				if (error > 1000.0 * tol_rel) {
					if (verbose) std::cout << "Residual too large, terminating early\n";
					return false;
				}
			}
		}
		return false;
	}
}

#endif // !GMRES_H
