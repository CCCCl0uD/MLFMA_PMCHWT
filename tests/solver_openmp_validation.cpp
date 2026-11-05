#include "solver_test_common.h"
#include <atomic>
#include <functional>

using namespace SolverTest;

static void validate_kernels(int n) {
    Vector a = known_solution(n), b(n), output(n), shadow(n), reference(n);
    long double conj_re = 0, conj_im = 0, bil_re = 0, bil_im = 0;
    double scale = 1;
    for (int i = 0; i < n; ++i) {
        b[i] = Complex((i % 2 ? -1 : 1) * (0.3 + i % 19), 0.1 * (i % 13));
        const Complex c = std::conj(a[i]) * b[i], d = a[i] * b[i];
        conj_re += c.real(); conj_im += c.imag();
        bil_re += d.real(); bil_im += d.imag();
        scale += std::abs(c);
    }
    const Complex conj_ref(static_cast<double>(conj_re), static_cast<double>(conj_im));
    const Complex bil_ref(static_cast<double>(bil_re), static_cast<double>(bil_im));
    require(std::abs(MGMRES::dot_conj(n, a.data(), b.data()) - conj_ref) < 1e-12 * scale,
        "conjugate dot mismatch");
    require(std::abs(MCGS::detail::dot_bilinear(n, a.data(), b.data()) - bil_ref) < 1e-12 * scale,
        "bilinear dot mismatch");
    require(std::abs(SolverVectorOps::norm2(n, a.data()) - norm(a)) < 1e-12 * (1 + norm(a)),
        "norm mismatch");
    SolverVectorOps::copy(n, a.data(), output.data());
    require(output == a, "copy mismatch");
    for (int i = 0; i < n; ++i) reference[i] = a[i] - b[i];
    double rn = SolverVectorOps::residual_norm(n, a.data(), b.data(), output.data(), shadow.data());
    require(output == reference && shadow == reference, "residual/shadow mismatch");
    require(std::abs(rn - norm(reference)) < 1e-12 * (1 + rn), "residual norm mismatch");
    const Complex alpha(0.7, -0.2);
    for (int i = 0; i < n; ++i) reference[i] -= alpha * b[i];
    rn = SolverVectorOps::update_residual_norm(n, output.data(), alpha, b.data());
    require(output == reference, "residual update mismatch");
    require(std::abs(rn - norm(reference)) < 1e-12 * (1 + rn), "updated norm mismatch");

    constexpr int m = 7;
    Vector v(static_cast<std::size_t>(n) * m), y(m);
    for (int j = 0; j < m; ++j) {
        y[j] = Complex(0.1 * j, -0.07);
        for (int i = 0; i < n; ++i) v[static_cast<std::size_t>(j) * n + i] = a[i] * Complex(j + 1, 0.3);
    }
    reference = output;
    for (int j = 0; j < m; ++j)
        for (int i = 0; i < n; ++i) reference[i] += y[j] * v[static_cast<std::size_t>(j) * n + i];
    SolverVectorOps::update_solution(n, m, output.data(), v.data(), y.data());
    require(relative_difference(output, reference) < 1e-14, "blocked solution update mismatch");
}

static void validate_solve(int n, int kind, bool nonzero_guess) {
    constexpr double tol = 1e-8;
    Matrix matrix;
    matrix.parallel_ax = true;
    Vector exact = known_solution(n), b(n), x(n);
    matrix(n, exact.data(), b.data());
    if (nonzero_guess)
        for (int i = 0; i < n; ++i) x[i] = Complex(0.12, -0.03);
    Vector reference = x;
    matrix.calls = 0;
    DiagonalPreconditioner preconditioner{n};
    auto ax = [&](int size, const Complex* input, Complex* output) { matrix(size, input, output); };
    bool success;
    if (kind == 0) success = MGMRES::mgmres(n, x.data(), b.data(), 100, 3, 0, tol, ax, false);
    else if (kind == 1) success = MCGS::CGS::solve(n, x.data(), b.data(), 100, tol, ax, false);
    else success = MCGS::CGS::solve(n, x.data(), b.data(), 100, tol, ax, preconditioner, false);
    const int calls = matrix.calls;
    require(success, "solver did not converge, kind=" + std::to_string(kind));
    require(!matrix.bad_context && !preconditioner.bad_context, "callback entered inside solver parallel region");
    require(kind != 2 || preconditioner.calls == calls - 1, "preconditioner call count mismatch");
    const double residual = true_residual(matrix, x, b);
    require(std::isfinite(residual) && residual < tol, "true residual exceeds tolerance");
    require(relative_difference(x, exact) < 10 * tol, "known solution mismatch");
#ifdef MLFMA_TEST_BASELINE
    Matrix baseline_matrix;
    baseline_matrix.parallel_ax = true;
    auto baseline_ax = [&](int size, const Complex* input, Complex* output) { baseline_matrix(size, input, output); };
    DiagonalPreconditioner baseline_preconditioner{n};
    bool baseline_success;
    if (kind == 0) baseline_success = SerialMGMRES::mgmres(n, reference.data(), b.data(), 100, 3, 0, tol, baseline_ax, false);
    else if (kind == 1) baseline_success = SerialMCGS::CGS::solve(n, reference.data(), b.data(), 100, tol, baseline_ax, false);
    else baseline_success = SerialMCGS::CGS::solve(n, reference.data(), b.data(), 100, tol, baseline_ax, baseline_preconditioner, false);
    require(success == baseline_success, "baseline status mismatch");
    require(calls == baseline_matrix.calls, "baseline callback count mismatch");
    require(relative_difference(x, reference) < 10 * tol, "baseline solution mismatch");
#endif
    std::cout << "solve n=" << n << " kind=" << kind << " guess=" << nonzero_guess
        << " calls=" << calls << " residual=" << residual << "\n";
}

static void validate_edges() {
    const int n = 32771;
    Vector zero(n), x(n);
    int calls = 0;
    auto identity = [&](int size, const Complex* input, Complex* output) {
        ++calls;
        std::copy_n(input, size, output);
    };
    require(MGMRES::mgmres(n, x.data(), zero.data(), 10, 3, 0, 1e-8, identity, false), "GMRES zero RHS");
    require(calls == 1, "GMRES initial convergence callback count");
    calls = 0;
    require(MCGS::CGS::solve(n, x.data(), zero.data(), 10, 1e-8, identity, false), "CGS zero RHS");
    require(calls == 1, "CGS initial convergence callback count");
    Vector exact = known_solution(n);
    x = exact;
    require(MGMRES::mgmres(n, x.data(), exact.data(), 10, 3, 0, 1e-8, identity, false), "GMRES exact initial guess");
    require(MCGS::CGS::solve(n, x.data(), exact.data(), 10, 1e-8, identity, false), "CGS exact initial guess");
    x.assign(n, 0);
    require(!MGMRES::mgmres(n, x.data(), exact.data(), 0, 3, 0, 1e-8, identity, false), "GMRES iteration limit");
    require(!MCGS::CGS::solve(n, x.data(), exact.data(), 0, 1e-8, identity, false), "CGS iteration limit");
    require(!MCGS::CGS::solve(0, x.data(), exact.data(), 2, 1e-8, identity, false), "CGS invalid dimension");
    // A nonzero isotropic shadow residual makes the bilinear rho exactly zero.
    Complex isotropic[2] = {Complex(1, 0), Complex(0, 1)}, small_x[2] = {};
    require(!MCGS::CGS::solve(2, small_x, isotropic, 3, 1e-8, identity, false), "CGS rho breakdown");
    auto zero_matrix = [](int size, const Complex*, Complex* out) { std::fill_n(out, size, Complex(0)); };
    require(!MCGS::CGS::solve(2, small_x, isotropic + 0, 3, 1e-8, zero_matrix, false), "CGS zero matrix breakdown");
    Complex rhs[2] = {Complex(1, 0), Complex(2, 0)};
    require(!MCGS::CGS::solve(2, small_x, rhs, 3, 1e-8, zero_matrix, false), "CGS alpha breakdown");
#ifdef MLFMA_TEST_BASELINE
    // Record existing happy-breakdown behavior rather than changing the algorithm.
    Complex baseline_x[2] = {};
    const bool current = MGMRES::mgmres(2, small_x, rhs, 2, 2, 0, 1e-8, identity, false);
    const bool baseline = SerialMGMRES::mgmres(2, baseline_x, rhs, 2, 2, 0, 1e-8, identity, false);
    require(current == baseline, "GMRES linear-dependence baseline status");
    for (int i = 0; i < 2; ++i) require(std::abs(small_x[i] - baseline_x[i]) < 1e-12, "GMRES linear-dependence baseline solution");
#endif
}

int main() {
    try {
#ifdef _OPENMP
        omp_set_dynamic(0);
        const std::vector<int> threads{1, 2, 4, 8, 12, 16};
#else
        const std::vector<int> threads{1};
#endif
        std::cout << std::scientific << std::setprecision(9);
        for (int count : threads) {
#ifdef _OPENMP
            omp_set_num_threads(count);
#endif
            std::cout << "threads=" << count << "\n";
            constexpr int largeSize = 32771;
#if defined(_OPENMP) && !defined(MLFMA_SOLVER_DISABLE_OPENMP)
            require(SolverVectorOps::use_parallel() == (count > 1), "parallel gate mismatch");
            std::atomic<int> workers{0};
            SolverVectorOps::for_each(largeSize, [&](int) {
                int observed = omp_get_num_threads();
                if (observed > workers.load()) workers.store(observed);
            });
            require(count == 1 || workers.load() > 1, "parallel kernel used only one worker");
#endif
            for (int repeat = 0; repeat < 3; ++repeat) {
                for (int n : {1, 7, largeSize}) validate_kernels(n);
                for (int n : {257, largeSize})
                    for (int kind = 0; kind < 3; ++kind)
                        for (bool guess : {false, true}) validate_solve(n, kind, guess);
                validate_edges();
            }
#ifdef _OPENMP
            // A single caller within an outer team must use serial vector kernels.
            bool nested_ok = false;
#pragma omp parallel num_threads(2)
            {
#pragma omp single
                {
                    Vector a = known_solution(largeSize), b(a.size());
                    nested_ok = !SolverVectorOps::use_parallel();
                    SolverVectorOps::copy(static_cast<int>(a.size()), a.data(), b.data());
                    nested_ok = nested_ok && a == b;
                }
            }
            require(nested_ok, "nested parallel fallback failed");
#endif
        }
        std::cout << "PASS: all solver and vector checks\n";
        return 0;
    }
    catch (const std::exception& e) {
        std::cerr << "FAIL: " << e.what() << "\n";
        return 1;
    }
}
