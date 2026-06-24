// CGS.h
#pragma once
#ifndef CGS_H
#define CGS_H

#include <algorithm>
#include <cmath>
#include <complex>
#include <iostream>
#include <type_traits>
#include <vector>

namespace MCGS {
	namespace detail {
		using Complex = std::complex<double>;

		// CGS.f90 中 rho 和 alpha 分母使用 sum(a*b)，不对 a 取共轭。
		inline Complex dot_bilinear(
			int n, const Complex* a, const Complex* b) {
			Complex sum(0.0, 0.0);
			for (int i = 0; i < n; ++i) {
				sum += a[i] * b[i];
			}
			return sum;
		}

		// 欧氏范数：对应 Fortran 的 sqrt(real(dot_product(a,a)))。
		inline double norm2(int n, const Complex* a) {
			double sum = 0.0;
			for (int i = 0; i < n; ++i) {
				sum += std::norm(a[i]);
			}
			return std::sqrt(sum);
		}

		struct NoPreconditioner {};

		template <typename Preconditioner>
		inline void apply_preconditioner(
			int n,
			Preconditioner& preconditioner,
			const Complex* input,
			Complex* output) {
			if constexpr (std::is_same_v<std::decay_t<Preconditioner>,
				NoPreconditioner>) {
				std::copy_n(input, n, output);
			}
			else {
				preconditioner.apply(input, output);
			}
		}
	} // namespace detail

	/**
	* 共轭梯度平方法（CGS）迭代求解器。
	*
	* AxFunc 的调用签名必须为：
	* void(int n, const std::complex<double>* x,
	*      std::complex<double>* w)
	*
	* 可选预处理子的接口必须为：
	* void apply(const std::complex<double>* p,
	*            std::complex<double>* P_)
	*/
	class CGS {
	public:
		/**
		* 无预处理版本。
		*
		* @return true 表示收敛；false 表示未收敛或发生算法击穿。
		*/
		template <typename AxFunc>
		static bool solve(
			int n,
			std::complex<double>* x,
			const std::complex<double>* b,
			int max_iter,
			double tol,
			AxFunc ax,
			bool verbose = true) {
			detail::NoPreconditioner no_preconditioner;
			return solve_impl(
				n, x, b, max_iter, tol, ax, no_preconditioner, verbose);
		}

		/**
		* 带预处理版本。preconditioner 可以保存矩阵维度等内部状态。
		*/
		template <typename AxFunc, typename Preconditioner>
		static bool solve(
			int n,
			std::complex<double>* x,
			const std::complex<double>* b,
			int max_iter,
			double tol,
			AxFunc ax,
			Preconditioner& preconditioner,
			bool verbose = true) {
			return solve_impl(
				n, x, b, max_iter, tol, ax, preconditioner, verbose);
		}

	private:
		template <typename AxFunc, typename Preconditioner>
		static bool solve_impl(
			int n,
			std::complex<double>* x,
			const std::complex<double>* b,
			int max_iter,
			double tol,
			AxFunc& ax,
			Preconditioner& preconditioner,
			bool verbose) {
			using Complex = std::complex<double>;
			constexpr double eps = 1.0e-16;

			if (n <= 0 || x == nullptr || b == nullptr ||
				max_iter < 0 || tol < 0.0) {
				return false;
			}

			// 工作向量命名与 CGS.f90 保持一致。
			std::vector<Complex> r(n), R_tilde(n), u(n), p(n);
			std::vector<Complex> P_(n), V_(n), q(n);
			std::vector<Complex> u_plus_q(n), U_(n), Q_(n);

			// 初始残差 r = b - A*x；影子残差 R~ = r0，此后不再改变。
			ax(n, x, V_.data());
			for (int i = 0; i < n; ++i) {
				r[i] = b[i] - V_[i];
				R_tilde[i] = r[i];
			}

			// 与 GMRES.h 一致：零右端时使用 1，避免相对残差除零。
			double bnorm = detail::norm2(n, b);
			if (bnorm < eps) {
				bnorm = 1.0;
			}

			double RSS = detail::norm2(n, r.data()) / bnorm;
			if (RSS < tol) {
				return true;
			}

			Complex rho_previous(0.0, 0.0);

			// 对应 CGS.f90 第 64–115 行。
			for (int iter = 1; iter <= max_iter; ++iter) {
				// rho(1) = sum(r * R_)
				const Complex rho =
					detail::dot_bilinear(n, r.data(), R_tilde.data());
				if (std::abs(rho) < eps) {
					if (verbose) {
						std::cout << "CGS breakdown: rho is zero at iteration "
							<< iter << "\n";
					}
					return false;
				}

				// 首次迭代：u = r, p = u；其后使用 beta 更新。
				if (iter == 1) {
					u = r;
					p = u;
				}
				else {
					const Complex beta = rho / rho_previous;
					for (int i = 0; i < n; ++i) {
						u[i] = r[i] + beta * q[i];
						p[i] = u[i] + beta * (q[i] + beta * p[i]);
					}
				}

				// P_ = M^{-1}p；无预处理时 P_ = p。
				detail::apply_preconditioner(
					n, preconditioner, p.data(), P_.data());

				// 第一次矩阵-向量乘：V_ = A*P_。
				ax(n, P_.data(), V_.data());

				const Complex alpha_denominator =
					detail::dot_bilinear(n, R_tilde.data(), V_.data());
				if (std::abs(alpha_denominator) < eps) {
					if (verbose) {
						std::cout
							<< "CGS breakdown: alpha denominator is zero at iteration "
							<< iter << "\n";
					}
					return false;
				}

				const Complex alpha = rho / alpha_denominator;
				for (int i = 0; i < n; ++i) {
					q[i] = u[i] - alpha * V_[i];
					u_plus_q[i] = u[i] + q[i];
				}

				// U_ = M^{-1}(u+q)；无预处理时 U_ = u+q。
				detail::apply_preconditioner(
					n, preconditioner, u_plus_q.data(), U_.data());

				for (int i = 0; i < n; ++i) {
					x[i] += alpha * U_[i];
				}

				// 第二次矩阵-向量乘：Q_ = A*U_。
				ax(n, U_.data(), Q_.data());
				for (int i = 0; i < n; ++i) {
					r[i] -= alpha * Q_[i];
				}

				// RSS = ||r|| / ||b||。
				RSS = detail::norm2(n, r.data()) / bnorm;
				if (verbose) {
					std::cout << "Iter " << iter << " Res: " << RSS << "\n";
				}

				if (RSS < tol) {
					if (verbose) {
						std::cout << "CGS converged! Final residual: "
							<< RSS << "\n";
					}
					return true;
				}

				// rho(2) = rho(1)，供下一次迭代计算 beta。
				rho_previous = rho;
			}

			if (verbose) {
				std::cout << "CGS reached max iterations, final residual: "
					<< RSS << "\n";
			}
			return false;
		}
	};
} // namespace MCGS

#endif // CGS_H