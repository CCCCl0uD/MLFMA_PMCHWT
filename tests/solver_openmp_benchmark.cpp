#include "solver_test_common.h"

using namespace SolverTest;
static volatile double checksum = 0;

template <typename Operation>
static double median_seconds(int repetitions, Operation operation) {
    operation(); // Warm up the runtime and data before recording samples.
    std::vector<double> samples;
    for (int repeat = 0; repeat < repetitions; ++repeat) {
        const auto start = Clock::now();
        operation();
        samples.push_back(std::chrono::duration<double>(Clock::now() - start).count());
    }
    std::sort(samples.begin(), samples.end());
    return samples[samples.size() / 2];
}

int main(int argc, char** argv) {
    try {
        const int n = argc > 1 ? std::stoi(argv[1]) : 100000;
        const int threads = argc > 2 ? std::stoi(argv[2]) : 4;
        const int repetitions = argc > 3 ? std::stoi(argv[3]) : 5;
        const int mr = argc > 4 ? std::stoi(argv[4]) : 50;
        require(n > 0 && threads > 0 && repetitions > 0 && mr > 0, "invalid benchmark parameters");
#ifdef _OPENMP
        omp_set_dynamic(0);
        omp_set_num_threads(threads);
#endif
        std::cout << std::setprecision(9);
        std::cout << "operation,n,mr,threads,mode,median_seconds,ax_seconds,ax_calls,true_residual\n";
        Vector exact = known_solution(n), b(n), x(n);
        Matrix matrix;
        matrix.parallel_ax = true;
        matrix(n, exact.data(), b.data());
        for (int kind = 0; kind < 2; ++kind) {
            const int modes =
#ifdef MLFMA_TEST_BASELINE
                2;
#else
                1;
#endif
            for (int mode = 0; mode < modes; ++mode) {
                std::vector<double> times, ax_times;
                int calls = 0;
                double residual = 0;
                for (int repeat = -1; repeat < repetitions; ++repeat) {
                    std::fill(x.begin(), x.end(), Complex(0));
                    matrix.calls = 0;
                    matrix.seconds = 0;
                    auto ax = [&](int size, const Complex* in, Complex* out) { matrix(size, in, out); };
                    const auto start = Clock::now();
                    bool success = false;
                    if (mode == 0) {
                        if (kind == 0) success = MGMRES::mgmres(n, x.data(), b.data(), 200, mr, 0, 1e-8, ax, false);
                        else success = MCGS::CGS::solve(n, x.data(), b.data(), 200, 1e-8, ax, false);
                    }
#ifdef MLFMA_TEST_BASELINE
                    else {
                        if (kind == 0) success = SerialMGMRES::mgmres(n, x.data(), b.data(), 200, mr, 0, 1e-8, ax, false);
                        else success = SerialMCGS::CGS::solve(n, x.data(), b.data(), 200, 1e-8, ax, false);
                    }
#endif
                    const double seconds = std::chrono::duration<double>(Clock::now() - start).count();
                    require(success && !matrix.bad_context, "benchmark solve failed");
                    calls = matrix.calls;
                    if (repeat >= 0) { times.push_back(seconds); ax_times.push_back(matrix.seconds); }
                    residual = true_residual(matrix, x, b);
                    require(residual < 1e-8, "benchmark true residual failed");
                    checksum += x[n / 2].real();
                }
                std::sort(times.begin(), times.end());
                std::sort(ax_times.begin(), ax_times.end());
                std::cout << (kind == 0 ? "GMRES" : "CGS") << ',' << n << ',' << mr << ',' << threads
                    << ',' << (mode == 0 ? "current" : "baseline") << ',' << times[times.size()/2]
                    << ',' << ax_times[ax_times.size()/2] << ',' << calls << ',' << residual << '\n';
            }
        }

        for (int mode = 0; mode < 2; ++mode) {
            const int batch = (std::max)(2, (std::min)(100, 1000000 / n));
            const double seconds = median_seconds(repetitions, [&] {
                for (int repeat = 0; repeat < batch; ++repeat) {
                    Complex value = 0;
                    if (mode == 0) value = SolverVectorOps::dot<true>(n, exact.data(), b.data());
                    else for (int i = 0; i < n; ++i) value += std::conj(exact[i]) * b[i];
                    checksum += value.real();
                }
            }) / batch;
            std::cout << "dot," << n << ',' << mr << ',' << threads << ','
                << (mode == 0 ? "current" : "baseline") << ',' << seconds << ",0,0,0\n";
        }

        // Exercise a full two-pass MGS column with mr preceding basis vectors.
        Vector basis(static_cast<std::size_t>(n) * mr), initial(n), w(n);
        for (int j = 0; j < mr; ++j) {
            const double factor = 1 / std::sqrt(static_cast<double>(n));
            for (int i = 0; i < n; ++i) basis[static_cast<std::size_t>(j) * n + i] =
                Complex(std::cos(i * (j + 1) * 0.031), std::sin(i * (j + 1) * 0.019)) * factor;
        }
        initial = exact;
        Vector current_result;
        for (int mode = 0; mode < 2; ++mode) {
            const double seconds = median_seconds(repetitions, [&] {
                w = initial;
                for (int pass = 0; pass < 2; ++pass) {
                    for (int j = 0; j < mr; ++j) {
                        const Complex* column = basis.data() + static_cast<std::size_t>(j) * n;
                        Complex value = 0;
                        if (mode == 0) {
                            value = SolverVectorOps::dot<true>(n, column, w.data());
                            SolverVectorOps::for_each(n, [&](int i) { w[i] -= value * column[i]; });
                        }
                        else {
                            for (int i = 0; i < n; ++i) value += std::conj(column[i]) * w[i];
                            for (int i = 0; i < n; ++i) w[i] -= value * column[i];
                        }
                    }
                }
                checksum += w[n / 2].real();
            });
            if (mode == 0) current_result = w;
            else require(relative_difference(current_result, w) < 1e-10, "MGS benchmark mismatch");
            std::cout << "MGS2," << n << ',' << mr << ',' << threads << ','
                << (mode == 0 ? "current" : "baseline") << ',' << seconds << ",0,0,0\n";
        }
        return 0;
    }
    catch (const std::exception& e) {
        std::cerr << e.what() << '\n';
        return 1;
    }
}
