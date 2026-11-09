#pragma once
#include <algorithm>
#include <cmath>
#include <complex>
#include <array>
#include <iterator>
#include <limits>
#include <map>
#include <set>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>
#include "AlgoMLFMM.h"
#include "Operator.h"
#include "ProcessFMM.h"
#include "ProcessMLFMM.h"
#include "RWGGenerator.h"

class CalderonGramLU {
public:
    using Cx = std::complex<double>;
    using Rows = std::vector<std::map<int, double>>;

private:
    Rows lu_;
    std::vector<int> pivots_;
    std::vector<double> scales_;

public:
    void factor(Rows rows)
    {
        lu_ = std::move(rows);
        const int n = static_cast<int>(lu_.size());
        pivots_.resize(n);
        scales_.resize(n);
        std::vector<std::set<int>> active(n);
        std::size_t inputNonzeros = 0;
        for (int r = 0; r < n; ++r) {
            double scale = 0.0;
            for (const auto& [col, value] : lu_[r])
                scale = std::max(scale, std::abs(value));
            if (!(scale > 0.0))
                throw std::runtime_error("Calderon: zero Gram row");
            scales_[r] = scale;
            auto& row = lu_[r];
            for (auto it = row.begin(); it != row.end(); ) {
                it->second /= scale;
                if (it->second == 0.0) {
                    it = row.erase(it);
                } else {
                    active[it->first].insert(r);
                    ++inputNonzeros;
                    ++it;
                }
            }
        }
        const double tol = 64.0 * std::numeric_limits<double>::epsilon();
        for (int k = 0; k < n; ++k) {
            int p = -1;
            double largest = 0.0;
            for (const int r : active[k]) {
                const double v = std::abs(lu_[r].find(k)->second);
                if (v > largest) {
                    p = r;
                    largest = v;
                }
            }
            if (!(largest > tol))
                throw std::runtime_error("Calderon: singular or ill-scaled Gram pivot");
            pivots_[k] = p;
            const auto& pivot = lu_[p];
            const auto diagonal = pivot.find(k);
            for (auto it = diagonal; it != pivot.end(); ++it)
                active[it->first].erase(p);
            for (const int r : active[k]) {
                auto& row = lu_[r];
                auto lower = row.find(k);
                const double multiplier = lower->second /= diagonal->second;
                for (auto it = std::next(diagonal); it != pivot.end(); ++it) {
                    const int j = it->first;
                    auto [entry, inserted] = row.try_emplace(j, 0.0);
                    entry->second -= multiplier * it->second;
                    if (entry->second == 0.0) {
                        row.erase(entry);
                        active[j].erase(r);
                    } else if (inserted) {
                        active[j].insert(r);
                    }
                }
            }
            active[k].clear();
        }
        std::size_t factorNonzeros = 0;
        for (const auto& row : lu_) factorNonzeros += row.size();
        std::cout << "Calderon Gram: " << inputNonzeros << " entries, "
                  << factorNonzeros << " LU entries\n";
    }

    void solve(const Cx* rhs, Cx* z) const
    {
        const int n = static_cast<int>(pivots_.size());
        for (int k = 0; k < n; ++k) {
            const int r = pivots_[k];
            Cx value = rhs[r] / scales_[r];
            const auto stop = lu_[r].lower_bound(k);
            for (auto it = lu_[r].begin(); it != stop; ++it)
                value -= it->second * z[it->first];
            z[k] = value;
        }
        for (int k = n - 1; k >= 0; --k) {
            const auto& row = lu_[pivots_[k]];
            const auto diagonal = row.find(k);
            for (auto it = std::next(diagonal); it != row.end(); ++it)
                z[k] -= it->second * z[it->first];
            z[k] /= diagonal->second;
        }
    }
};

class CalderonGeometry {
public:
    using Cx = CalderonGramLU::Cx;
    using Rows = CalderonGramLU::Rows;

    std::vector<Point> finePoints;
    std::vector<FaceElement> fineFaces;
    std::vector<RWGBase> fineRwgs;
    Rows C;
    CalderonGramLU gramLU;
    const int n;

    CalderonGeometry(const CalderonGeometry&) = delete;
    CalderonGeometry& operator=(const CalderonGeometry&) = delete;

    CalderonGeometry(const std::vector<Point>& points,
                     const std::vector<FaceElement>& faces,
                     const std::vector<RWGBase>& rwgs,
                     const RWGGenerationReport& coarseReport)
        : n(static_cast<int>(rwgs.size()))
    {
        if (coarseReport.components.size() != 1)
            throw std::invalid_argument("Calderon v1 expects one closed interface");

        using Edge = std::pair<int, int>;
        const auto key = [](int a, int b) -> Edge {
            return {std::min(a, b), std::max(a, b)};
        };
        const auto crossPoint = [](const Point& a, const Point& b) -> Point {
            return {0, a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z,
                       a.x * b.y - a.y * b.x};
        };
        const auto dotPoint = [](const Point& a, const Point& b) {
            return a.x * b.x + a.y * b.y + a.z * b.z;
        };

        const int nv = static_cast<int>(points.size());
        const int nt = static_cast<int>(faces.size());
        const Point& origin = *faces.front().vertex1;
        double volume6 = 0.0;
        for (const auto& f : faces)
            volume6 += dotPoint(*f.vertex1 - origin,
                crossPoint(*f.vertex2 - origin, *f.vertex3 - origin));
        if (!(volume6 > 0.0))
            throw std::invalid_argument("Calderon: outward winding required");

        finePoints.resize(points.size() + rwgs.size() + faces.size());
        for (int i = 0; i < nv; ++i) finePoints[i] = points[i];
        std::map<Edge, int> midpoint;
        for (const auto& e : rwgs) {
            const int a = static_cast<int>(e.commonEdgeVertex1 - points.data());
            const int b = static_cast<int>(e.commonEdgeVertex2 - points.data());
            const int m = nv + e.rwgid;
            finePoints[m] = (points[a] + points[b]) / 2.0;
            midpoint.emplace(key(a, b), m);
        }
        for (int t = 0; t < nt; ++t)
            finePoints[nv + n + t] = faces[t].triangleCenter;
        for (std::size_t i = 0; i < finePoints.size(); ++i)
            finePoints[i].id = static_cast<int>(i) + 1;

        fineFaces.reserve(6 * faces.size());
        std::vector<int> owner;
        owner.reserve(6 * faces.size());
        std::vector<int> cellSize(nv, 0);
        const auto addFace = [&](int v, int a, int b) {
            FaceElement f{};
            f.faceId = static_cast<int>(fineFaces.size()) + 1;
            f.vertex1 = &finePoints[v];
            f.vertex2 = &finePoints[a];
            f.vertex3 = &finePoints[b];
            const Point normal = crossPoint(*f.vertex2 - *f.vertex1,
                                            *f.vertex3 - *f.vertex1);
            const double area2 = std::hypot(normal.x, normal.y, normal.z);
            f.triangularArea = area2 / 2.0;
            f.externalNormalVector = normal / area2;
            f.triangleCenter = (*f.vertex1 + *f.vertex2 + *f.vertex3) / 3.0;
            fineFaces.push_back(f);
            owner.push_back(v);
            ++cellSize[v];
        };
        for (int t = 0; t < nt; ++t) {
            const auto& f = faces[t];
            const int vertices[] = {
                static_cast<int>(f.vertex1 - points.data()),
                static_cast<int>(f.vertex2 - points.data()),
                static_cast<int>(f.vertex3 - points.data())
            };
            const int center = nv + n + t;
            for (int k = 0; k < 3; ++k) {
                const int v = vertices[k];
                addFace(v, midpoint.at(key(v, vertices[(k + 1) % 3])), center);
                addFace(v, center, midpoint.at(key(v, vertices[(k + 2) % 3])));
            }
        }

        RWGGenerationReport report;
        if (!RWGGenerator::GenerateRWG(fineFaces, finePoints, fineRwgs, RWGOptions{}, report))
            throw std::runtime_error("Calderon: refined RWG generation failed: " + report.error);

        std::map<Edge, int> fineEdge;
        std::vector<std::array<int, 2>> adjacent(fineRwgs.size());
        std::vector<std::vector<std::pair<int, int>>> inside(fineFaces.size());
        for (const auto& f : fineRwgs) {
            const int e = f.rwgid;
            fineEdge.emplace(key(f.commonEdgeVertex1->id - 1,
                                 f.commonEdgeVertex2->id - 1), e);
            const int t = f.positiveFace->faceId - 1;
            const int u = f.negativeFace->faceId - 1;
            adjacent[e] = {t, u};
            if (owner[t] == owner[u]) {
                inside[t].emplace_back(u, e);
                inside[u].emplace_back(t, e);
            }
        }

        C.resize(fineRwgs.size());
        for (const auto& e : rwgs) {
            const int col = e.rwgid;
            const int a = static_cast<int>(e.commonEdgeVertex1 - points.data());
            const int b = static_cast<int>(e.commonEdgeVertex2 - points.data());
            const int mid = nv + col;
            std::map<int, double> fixedOut;
            const auto writeFlux = [&](int edge, double q) {
                if (q != 0.0)
                    C[edge].emplace(col, q * e.edgeLength / fineRwgs[edge].edgeLength);
            };
            for (const auto* parent : {e.positiveFace, e.negativeFace}) {
                const int center = nv + n + static_cast<int>(parent - faces.data());
                const int edge = fineEdge.at(key(mid, center));
                const auto [t, u] = adjacent[edge];
                const double q = owner[t] == a ? 0.5 : -0.5;
                writeFlux(edge, q);
                fixedOut[t] += q;
                fixedOut[u] -= q;
            }
            for (const int v : {a, b}) {
                const int cut = fineEdge.at(key(v, mid));
                int t = adjacent[cut][0];
                int previous = cut;
                const double target = (v == a ? 1.0 : -1.0) / cellSize[v];
                double incoming = 0.0;
                // Each refined face has two edges in its cell. The cut opens the ring.
                for (int step = 1; step < cellSize[v]; ++step) {
                    const auto& edges = inside[t];
                    const auto [next, edge] = edges.front().second == previous ? edges.back() : edges.front();
                    const auto known = fixedOut.find(t);
                    const double outward = known == fixedOut.end() ? 0.0 : known->second;
                    const double q = target + incoming - outward;
                    writeFlux(edge, adjacent[edge][0] == t ? q : -q);
                    incoming = q;
                    t = next;
                    previous = edge;
                }
            }
        }

        struct LocalBasis {
            int id;
            const Point* free;
            double scale;
        };
        const auto faceBasis = [](const std::vector<FaceElement>& mesh,
                                  const std::vector<RWGBase>& bases) {
            std::vector<std::vector<LocalBasis>> local(mesh.size());
            for (const auto& e : bases) {
                local[e.positiveFace - mesh.data()].push_back(
                    {e.rwgid, e.freeVertexPositive,
                     e.edgeLength / (2.0 * e.positiveFace->triangularArea)});
                local[e.negativeFace - mesh.data()].push_back(
                    {e.rwgid, e.freeVertexNegative,
                     -e.edgeLength / (2.0 * e.negativeFace->triangularArea)});
            }
            return local;
        };
        const auto coarseBasis = faceBasis(faces, rwgs);
        const auto fineBasis = faceBasis(fineFaces, fineRwgs);
        Rows G(n);
        for (std::size_t t = 0; t < fineFaces.size(); ++t) {
            const auto& f = fineFaces[t];
            for (const auto& i : coarseBasis[t / 6]) {
                const Point rotated = crossPoint(f.externalNormalVector,
                                                (f.triangleCenter - *i.free) * i.scale);
                for (const auto& a : fineBasis[t]) {
                    const double mass = f.triangularArea * a.scale *
                        dotPoint(rotated, f.triangleCenter - *a.free);
                    for (const auto& [j, coeff] : C[a.id])
                        G[i.id][j] += mass * coeff;
                }
            }
        }
        gramLU.factor(std::move(G));
    }
};

class CalderonV1 {
    using Cx = CalderonGeometry::Cx;
    using Vec = std::vector<Cx>;
    using Angles = std::vector<std::vector<double>>;
    using Waves = std::vector<std::vector<kp_Point>>;
    using Work = std::vector<std::vector<std::vector<Complex3D>>>;

    const CalderonGeometry& geo_;
    const Cx kappa_;
    OCTree tree_;
    int rows_, span_;
    Angles phi_, weights_;
    std::vector<double> phiWeights_;
    Waves kp_;
    std::vector<std::vector<InterpData>> interp_;
    std::vector<std::vector<std::vector<Cx>>> transfer_;
    std::vector<std::vector<Complex3D>> source_, test_;
    Work outgoing_, incoming_;
    std::vector<std::vector<std::pair<int, Cx>>> near_;
    Vec bc_, fineX_, fineY_;

    void applyFine(const Cx* x, Cx* y)
    {
        const int depth = tree_.maxLevel_;
        const auto& levels = tree_.octreeNodes_;
        const auto& leaves = levels[depth];
        ProcessMLFMM::clearSmBm(outgoing_, incoming_);
        ProcessMLFMM::mexp(outgoing_, leaves, source_, kp_[0], x, 0, depth);
        ProcessMLFMM::m2m(outgoing_, levels, kp_, interp_, depth, span_, kappa_);
        ProcessMLFMM::m2l(outgoing_, incoming_, levels, transfer_, kp_,
                         weights_, phiWeights_, phi_, depth, span_);
        ProcessMLFMM::l2l(incoming_, levels, kp_, phi_, interp_, depth, span_, kappa_);
        const Cx scale = kappa_ * kappa_ / (16.0 * Pi * Pi);
        const int leafCount = static_cast<int>(leaves.size());
#pragma omp parallel for schedule(dynamic)
        for (int t = 0; t < leafCount; ++t) {
            const auto& local = incoming_[depth - 2][t];
            for (const auto* f : leaves[t]->rwgIndices) {
                const int r = f->rwgid;
                Cx far{};
                for (std::size_t k = 0; k < local.size(); ++k)
                    far += dot(test_[r][k], local[k]);
                Cx value = scale * far;
                for (const auto& entry : near_[r])
                    value += entry.second * x[entry.first];
                y[r] = value;
            }
        }
    }

public:
    CalderonV1(const CalderonGeometry& geo, const EMSource& wave,
               const gaussPoints& gp, const std::string& nas)
        : geo_(geo), kappa_(wave.k1() * Cx(1.0, -0.01)),
          tree_(geo.fineRwgs, wave.wavelength(), nas),
          bc_(geo.n), fineX_(geo.fineRwgs.size()), fineY_(geo.fineRwgs.size())
    {
        const int depth = tree_.maxLevel_;
        const auto& levels = tree_.octreeNodes_;
        std::vector<int> orders, unusedOrders;
        Angles theta, unusedWeights, unusedTheta, unusedPhi;
        std::vector<double> unusedPhiWeights;
        Waves unusedKp;
        AlgoMLFMM::computeBase_MLFMM(
            levels, std::abs(kappa_), 0.0, 0, depth,
            static_cast<int>(geo.fineRwgs.size()), orders, unusedOrders, rows_, span_,
            weights_, unusedWeights, phiWeights_, unusedPhiWeights,
            theta, unusedTheta, phi_, unusedPhi, kp_, unusedKp);
        AlgoMLFMM::computeTransferFactorMatrix(
            transfer_, tree_.octreeNodesDRvec_, kp_, orders, depth, kappa_);
        AlgoMLFMM::computeInterpolationMatrix(interp_, kp_, theta, phi_, depth);
        const auto& leaves = levels[depth];
        ProcessFMM::computeVsmi_Vfmj_Die(
            source_, test_, leaves, kp_[0], kappa_, rows_, gp, 1.0, false);

        outgoing_.resize(span_);
        incoming_.resize(span_);
        for (int level = 2; level <= depth; ++level) {
            const auto nodes = levels[level].size();
            const auto samples = kp_[depth - level].size();
            outgoing_[level - 2].assign(nodes, std::vector<Complex3D>(samples));
            incoming_[level - 2].assign(nodes, std::vector<Complex3D>(samples));
        }
        near_.resize(rows_);
        const int leafCount = static_cast<int>(leaves.size());
#pragma omp parallel for schedule(dynamic)
        for (int t = 0; t < leafCount; ++t) {
            const auto* target = leaves[t];
            for (const auto* f : target->rwgIndices) {
                auto& entries = near_[f->rwgid];
                const auto fill = [&](const OCTree::Node* source) {
                    for (const auto* g : source->rwgIndices) {
                        Cx value{};
                        OLK::L_operator(f, g, kappa_, gp, value);
                        entries.emplace_back(g->rwgid, value);
                    }
                };
                fill(target);
                for (const auto* s : target->nearNeighbors) fill(s);
            }
        }
    }

    void apply(const Cx* input, Cx* output)
    {
        const int n = geo_.n;
        for (int channel = 0; channel < 2; ++channel) {
            const auto* r = input + channel * n;
            auto* q = output + channel * n;
            geo_.gramLU.solve(r, bc_.data());
            std::fill(fineX_.begin(), fineX_.end(), Cx{});
            for (int a = 0; a < rows_; ++a)
                for (const auto& [j, coeff] : geo_.C[a])
                    fineX_[a] += coeff * bc_[j];
            applyFine(fineX_.data(), fineY_.data());
            std::fill(q, q + n, Cx{});
            for (int a = 0; a < rows_; ++a)
                for (const auto& [j, coeff] : geo_.C[a])
                    q[j] += coeff * fineY_[a];
        }
    }
};
