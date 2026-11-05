#pragma once

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstddef>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace SolverVectorOps {
    using Complex = std::complex<double>;

    inline bool use_parallel() {
#if defined(_OPENMP) && !defined(MLFMA_SOLVER_DISABLE_OPENMP)
        return omp_get_max_threads() > 1 && !omp_in_parallel();
#else
        return false;
#endif
    }

    // Each invocation owns one output element; callbacks must not mutate shared
    // state. Explicit serial branches avoid serialized OpenMP region overhead.
    template <typename Function>
    inline void for_each(int n, const Function& operation) {
#if defined(_OPENMP) && !defined(MLFMA_SOLVER_DISABLE_OPENMP)
        if (use_parallel()) {
#pragma omp parallel for schedule(static)
            for (int i = 0; i < n; ++i) operation(i);
            return;
        }
#endif
        for (int i = 0; i < n; ++i) operation(i);
    }

    template <bool Conjugate>
    inline Complex dot(int n, const Complex* a, const Complex* b) {
#if defined(_OPENMP) && !defined(MLFMA_SOLVER_DISABLE_OPENMP)
        if (use_parallel()) {
            double real_sum = 0.0, imag_sum = 0.0;
#pragma omp parallel for schedule(static) reduction(+:real_sum,imag_sum)
            for (int i = 0; i < n; ++i) {
                const Complex term = (Conjugate ? std::conj(a[i]) : a[i]) * b[i];
                real_sum += term.real();
                imag_sum += term.imag();
            }
            return Complex(real_sum, imag_sum);
        }
#endif
        Complex sum(0.0, 0.0);
        for (int i = 0; i < n; ++i)
            sum += (Conjugate ? std::conj(a[i]) : a[i]) * b[i];
        return sum;
    }

    inline double norm2(int n, const Complex* a) {
        double sum = 0.0;
#if defined(_OPENMP) && !defined(MLFMA_SOLVER_DISABLE_OPENMP)
        if (use_parallel()) {
#pragma omp parallel for schedule(static) reduction(+:sum)
            for (int i = 0; i < n; ++i) sum += std::norm(a[i]);
            return std::sqrt(sum);
        }
#endif
        for (int i = 0; i < n; ++i) sum += std::norm(a[i]);
        return std::sqrt(sum);
    }

    inline void copy(int n, const Complex* input, Complex* output) {
        for_each(n, [=](int i) { output[i] = input[i]; });
    }

    // Residual and optional shadow residual are complete before returning.
    inline double residual_norm(int n, const Complex* b, const Complex* ax,
        Complex* residual, Complex* shadow = nullptr) {
        double sum = 0.0;
#if defined(_OPENMP) && !defined(MLFMA_SOLVER_DISABLE_OPENMP)
        if (use_parallel()) {
#pragma omp parallel for schedule(static) reduction(+:sum)
            for (int i = 0; i < n; ++i) {
                residual[i] = b[i] - ax[i];
                if (shadow) shadow[i] = residual[i];
                sum += std::norm(residual[i]);
            }
            return std::sqrt(sum);
        }
#endif
        for (int i = 0; i < n; ++i) {
            residual[i] = b[i] - ax[i];
            if (shadow) shadow[i] = residual[i];
            sum += std::norm(residual[i]);
        }
        return std::sqrt(sum);
    }

    inline double update_residual_norm(int n, Complex* residual,
        Complex alpha, const Complex* q) {
        double sum = 0.0;
#if defined(_OPENMP) && !defined(MLFMA_SOLVER_DISABLE_OPENMP)
        if (use_parallel()) {
#pragma omp parallel for schedule(static) reduction(+:sum)
            for (int i = 0; i < n; ++i) {
                residual[i] -= alpha * q[i];
                sum += std::norm(residual[i]);
            }
            return std::sqrt(sum);
        }
#endif
        for (int i = 0; i < n; ++i) {
            residual[i] -= alpha * q[i];
            sum += std::norm(residual[i]);
        }
        return std::sqrt(sum);
    }

    inline void update_solution(int n, int m, Complex* x,
        const Complex* v, const Complex* y) {
#if defined(_OPENMP) && !defined(MLFMA_SOLVER_DISABLE_OPENMP)
        if (use_parallel()) {
            constexpr int block_size = 256;
            const int blocks = n / block_size + (n % block_size != 0);
#pragma omp parallel for schedule(static)
            for (int block = 0; block < blocks; ++block) {
                const int begin = block * block_size;
                const int end = begin + (std::min)(block_size, n - begin);
                for (int j = 0; j < m; ++j) {
                    const Complex* column = v + static_cast<std::size_t>(j) * n;
                    for (int i = begin; i < end; ++i) x[i] += y[j] * column[i];
                }
            }
            return;
        }
#endif
        // Preserve the original accumulation and traversal order on fallback.
        for (int j = 0; j < m; ++j) {
            const Complex* column = v + static_cast<std::size_t>(j) * n;
            for (int i = 0; i < n; ++i) x[i] += y[j] * column[i];
        }
    }
} // namespace SolverVectorOps
