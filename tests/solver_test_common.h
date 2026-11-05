#pragma once

#include <algorithm>
#include <chrono>
#include <cmath>
#include <complex>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>
#include "../Code/GMRES.h"
#include "../Code/CGS.h"

// The runner extracts these byte-for-byte from the selected Git reference.
// Renaming namespaces lets both implementations run in the same process.
#ifdef MLFMA_TEST_BASELINE
#undef GMRES_H
#define MGMRES SerialMGMRES
#include "../output/solver_openmp_baseline/GMRES.h"
#undef MGMRES
#undef CGS_H
#define MCGS SerialMCGS
#include "../output/solver_openmp_baseline/CGS.h"
#undef MCGS
#endif

namespace SolverTest {
    using Complex = std::complex<double>;
    using Vector = std::vector<Complex>;
    using Clock = std::chrono::steady_clock;

    inline void require(bool condition, const std::string& message) {
        if (!condition) throw std::runtime_error(message);
    }

    inline double norm(const Vector& x) {
        long double sum = 0;
        for (Complex z : x) sum += static_cast<long double>(z.real()) * z.real()
            + static_cast<long double>(z.imag()) * z.imag();
        return std::sqrt(static_cast<double>(sum));
    }

    inline double relative_difference(const Vector& a, const Vector& b) {
        Vector difference(a.size());
        for (std::size_t i = 0; i < a.size(); ++i) difference[i] = a[i] - b[i];
        return norm(difference) / (std::max)(norm(b), 1e-30);
    }

    inline Vector known_solution(int n) {
        Vector result(n);
        for (int i = 0; i < n; ++i)
            result[i] = Complex(0.8 + 0.2 * std::sin(i * 0.17),
                0.15 * std::cos(i * 0.13));
        return result;
    }

    // Non-Hermitian, strictly diagonally dominant complex tridiagonal matrix.
    // This callback has a fixed threading policy independent of solver kernels.
    struct Matrix {
        int calls = 0;
        double seconds = 0;
        bool bad_context = false;
        bool parallel_ax = false;

        static Complex diagonal(int i) {
            return Complex(2.0 + (i % 97) / 97.0, 0.08 + (i % 11) * 0.003);
        }

        static Complex element(int i, int n, const Complex* x) {
            Complex value = diagonal(i) * x[i];
            if (i > 0) value += Complex(-0.22, 0.04) * x[i - 1];
            if (i + 1 < n) value += Complex(0.11, -0.03) * x[i + 1];
            return value;
        }

        void operator()(int n, const Complex* x, Complex* y) {
            const auto start = Clock::now();
#ifdef _OPENMP
            bad_context = bad_context || omp_in_parallel() != 0;
            if (parallel_ax) {
#pragma omp parallel for schedule(static)
                for (int i = 0; i < n; ++i) y[i] = element(i, n, x);
            }
            else
#endif
            {
                for (int i = 0; i < n; ++i) y[i] = element(i, n, x);
            }
            ++calls;
            seconds += std::chrono::duration<double>(Clock::now() - start).count();
        }
    };

    struct DiagonalPreconditioner {
        int n;
        int calls = 0;
        bool bad_context = false;
        void apply(const Complex* x, Complex* y) {
#ifdef _OPENMP
            bad_context = bad_context || omp_in_parallel() != 0;
#endif
            ++calls;
            for (int i = 0; i < n; ++i) y[i] = x[i] / Matrix::diagonal(i);
        }
    };

    inline double true_residual(Matrix& matrix, const Vector& x, const Vector& b) {
        Vector r(x.size());
        matrix(static_cast<int>(x.size()), x.data(), r.data());
        for (std::size_t i = 0; i < r.size(); ++i) r[i] = b[i] - r[i];
        const double bnorm = norm(b);
        return norm(r) / (bnorm < 1e-16 ? 1 : bnorm);
    }
}
