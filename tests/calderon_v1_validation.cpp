#include "../Code/CalderonV1.h"
#include "../Code/GMRES.h"
#include "../Code/NASReader.h"
#include "../Code/Zmartix.h"

#include <random>

int main(int argc, char* argv[]) {
    if (argc != 3) {
        std::cerr << "Usage: calderon_v1_validation <closed sphere NAS> <frequency Hz>\n";
        return 1;
    }
    using Complex = std::complex<double>;
    using Vector = std::vector<Complex>;
    omp_set_num_threads(2);

    std::vector<Point> points;
    std::vector<FaceElement> faces;
    std::vector<RWGBase> rwgs;
    MeshReadReport meshReport;
    RWGGenerationReport rwgReport;
    if (readNasData(argv[1], points, faces, meshReport) != 0 ||
        !RWGGenerator::GenerateRWG(faces, points, rwgs, RWGOptions{true}, rwgReport))
        return 1;

    CalderonGeometry geometry(points, faces, rwgs, rwgReport);
    const int n = geometry.n;
    const int nf = static_cast<int>(geometry.fineRwgs.size());
    std::vector<int> cellSize(points.size());
    for (const auto& face : geometry.fineFaces) ++cellSize[face.vertex1->id - 1];
    for (const auto& coarse : rwgs) {
        std::vector<double> divergence(geometry.fineFaces.size());
        for (const auto& fine : geometry.fineRwgs) {
            const auto& row = geometry.C[fine.rwgid];
            const auto coefficient = row.find(coarse.rwgid);
            if (coefficient == row.end()) continue;
            const double flux = coefficient->second * fine.edgeLength / coarse.edgeLength;
            divergence[fine.positiveFace->faceId - 1] += flux;
            divergence[fine.negativeFace->faceId - 1] -= flux;
        }
        const int a = static_cast<int>(coarse.commonEdgeVertex1 - points.data());
        const int b = static_cast<int>(coarse.commonEdgeVertex2 - points.data());
        for (std::size_t t = 0; t < divergence.size(); ++t) {
            const int owner = geometry.fineFaces[t].vertex1->id - 1;
            const double expected = owner == a ? 1.0 / cellSize[a]
                : owner == b ? -1.0 / cellSize[b] : 0.0;
            if (std::abs(divergence[t] - expected) > 1e-12)
                throw std::runtime_error("BC divergence does not match its dual cells");
        }
    }

    EMSource wave(std::stod(argv[2]), 0.0);
    wave.initDielectric({3.0, -0.1}, {1.0, 0.0});
    const gaussPoints gp(7);
    const Complex kappa = wave.k1() * Complex(1.0, -0.01);
    const auto& fine = geometry.fineRwgs;
    std::vector<Vector> denseL(nf, Vector(nf));
#pragma omp parallel for schedule(dynamic)
    for (int i = 0; i < nf; ++i)
        for (int j = 0; j < nf; ++j)
            OLK::L_operator(&fine[i], &fine[j], kappa, gp, denseL[i][j]);

    const auto denseQ = [&](const Complex* input, Complex* output) {
        Vector z(n), lifted(nf);
        std::fill_n(output, 2 * n, Complex{});
        for (int channel = 0; channel < 2; ++channel) {
            geometry.gramLU.solve(input + channel * n, z.data());
            for (int i = 0; i < nf; ++i) {
                lifted[i] = {};
                for (const auto& [j, value] : geometry.C[i]) lifted[i] += value * z[j];
            }
            for (int i = 0; i < nf; ++i) {
                Complex value{};
                for (int j = 0; j < nf; ++j) value += denseL[i][j] * lifted[j];
                for (const auto& [j, coefficient] : geometry.C[i])
                    output[channel * n + j] += coefficient * value;
            }
        }
    };

    std::mt19937 random(23);
    std::uniform_real_distribution<double> uniform(-1.0, 1.0);
    const int size = 2 * n;
    Vector expected(size), rhs(size), leftRhs(size), fast(size), raw(size), solution(size);
    for (auto& value : expected) value = {uniform(random), uniform(random)};
    Vector gramSolution(n), gramProduct(n);
    geometry.gramLU.solve(expected.data(), gramSolution.data());
    for (const auto& edge : fine) {
        Complex lifted{};
        for (const auto& [j, coefficient] : geometry.C[edge.rwgid])
            lifted += coefficient * gramSolution[j];
        for (const bool positive : {true, false}) {
            const auto* face = positive ? edge.positiveFace : edge.negativeFace;
            const auto* free = positive ? edge.freeVertexPositive : edge.freeVertexNegative;
            const auto* parent = &faces[(face->faceId - 1) / 6];
            const Point& normal = face->externalNormalVector;
            const double direction[3] = {normal.x, normal.y, normal.z};
            for (const auto& coarse : rwgs) {
                const bool plus = coarse.positiveFace == parent;
                if (!plus && coarse.negativeFace != parent) continue;
                const auto* coarseFree = plus ? coarse.freeVertexPositive : coarse.freeVertexNegative;
                const double scale = (positive == plus ? 1.0 : -1.0) * edge.edgeLength
                    * coarse.edgeLength / (12.0 * parent->triangularArea);
                for (const auto* vertex : {face->vertex1, face->vertex2, face->vertex3}) {
                    const Point sample = (face->triangleCenter + *vertex) / 2.0;
                    const Point coarseValue = sample - *coarseFree;
                    const Point fineValue = sample - *free;
                    const double test[3] = {coarseValue.x, coarseValue.y, coarseValue.z};
                    const double trial[3] = {fineValue.x, fineValue.y, fineValue.z};
                    gramProduct[coarse.rwgid] += lifted * scale * my_math::cross_dot(direction, test, trial);
                }
            }
        }
    }
    const double gramError = SolverVectorOps::residual_norm(n, expected.data(), gramProduct.data(), raw.data())
        / SolverVectorOps::norm2(n, expected.data());
    if (!(gramError <= 1e-12))
        throw std::runtime_error("Gram solve disagrees with three-point quadrature");

    std::vector<Vector> denseA;
    Zmartix::computeZ_die_pmchwt(rwgs, n, wave, gp, denseA, false);
    const auto physicalProduct = [&](const Complex* input, Complex* output) {
        for (int i = 0; i < size; ++i) {
            output[i] = {};
            for (int j = 0; j < size; ++j) output[i] += denseA[i][j] * input[j];
        }
    };
    physicalProduct(expected.data(), rhs.data());
    denseQ(rhs.data(), leftRhs.data());
    const auto leftProduct = [&](int, const Complex* input, Complex* output) {
        physicalProduct(input, raw.data());
        denseQ(raw.data(), output);
    };
    if (!MGMRES::mgmres(size, solution.data(), leftRhs.data(), size, size,
            0.0, 1e-8, leftProduct, false))
        throw std::runtime_error("Dense Calderon system did not converge");
    physicalProduct(solution.data(), raw.data());
    const double residual = SolverVectorOps::residual_norm(size, rhs.data(), raw.data(), raw.data())
        / SolverVectorOps::norm2(size, rhs.data());
    const double solutionError = SolverVectorOps::residual_norm(size, expected.data(), solution.data(), raw.data())
        / SolverVectorOps::norm2(size, expected.data());

    CalderonV1 fastQ(geometry, wave, gp, argv[1]);
    fastQ.apply(rhs.data(), fast.data());
    const double productError = SolverVectorOps::residual_norm(size, leftRhs.data(), fast.data(), raw.data())
        / SolverVectorOps::norm2(size, leftRhs.data());
    std::cout << std::scientific << "Gram relative residual: " << gramError
        << "\nDense original relative residual: " << residual
        << "\nDense current relative error: " << solutionError
        << "\nFast Q relative error: " << productError << '\n';
    return residual <= 1e-6 && solutionError <= 1e-6 && productError <= 1e-5 ? 0 : 1;
}
