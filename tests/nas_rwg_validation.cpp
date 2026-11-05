#include "../Code/OCTree.h"
#include "../Code/RWGGenerator.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <set>
#include <sstream>
#include <string>
#include <vector>

namespace fs = std::filesystem;

void require(bool condition, const std::string& message) {
    if (!condition) {
        std::cerr << "FAIL: " << message << '\n';
        std::exit(1);
    }
}

struct Mesh {
    std::vector<Point> points;
    std::vector<FaceElement> faces;
    MeshReadReport source;
    std::vector<RWGBase> bases;
};

std::string tetraPoints(int base, double x) {
    std::ostringstream out;
    out << "GRID," << base + 1 << ",," << x << ",0,0\n"
        << "GRID," << base + 2 << ",," << x + 1 << ",0,0\n"
        << "GRID," << base + 3 << ",," << x << ",1,0\n"
        << "GRID," << base + 4 << ",," << x << ",0,1\n";
    return out.str();
}

std::string tetraFaces(int base, int eid) {
    std::ostringstream out;
    const int indices[4][3] = {{1, 3, 2}, {1, 2, 4}, {2, 3, 4}, {3, 1, 4}};
    for (int i = 0; i < 4; ++i) {
        // One closed object deliberately spans four different PIDs.
        out << "CTRIA3," << eid + i << ',' << 100 + i;
        for (int v : indices[i]) out << ',' << base + v;
        out << '\n';
    }
    return out.str();
}

std::string largePoint(int id, double x, double y, double z, const std::string& label) {
    std::ostringstream out;
    out << std::scientific << std::setprecision(8)
        << std::left << std::setw(8) << "GRID*" << std::right
        << std::setw(16) << id << std::setw(16) << 0
        << std::setw(16) << x << std::setw(16) << y << std::setw(8) << label
        << "\n$ comment between the two halves\n\n"
        << std::left << std::setw(8) << ("*" + label) << std::right
        << std::setw(16) << z << '\n';
    return out.str();
}

fs::path writeFixture(const fs::path& dir, const std::string& name, const std::string& text) {
    const auto path = dir / (name + ".nas");
    std::ofstream out(path);
    out << text;
    out.close();
    require(bool(out), "fixture write: " + name);
    return path;
}

void readMesh(Mesh& mesh, const fs::path& path) {
    require(readNasData(path.string(), mesh.points, mesh.faces, mesh.source) == 0,
        "read: " + path.string());
    std::set<const Point*> owned;
    for (const auto& p : mesh.points) owned.insert(&p);
    for (const auto& face : mesh.faces) {
        require(owned.count(face.vertex1) && owned.count(face.vertex2) && owned.count(face.vertex3),
            "face pointers belong to final points storage");
        require(face.triangularArea > 0 && std::isfinite(face.triangularArea), "valid face area");
        const auto center = CalculateTriangleCenter(*face.vertex1, *face.vertex2, *face.vertex3);
        require(center.x == face.triangleCenter.x && center.y == face.triangleCenter.y &&
            center.z == face.triangleCenter.z, "center uses the referenced final vertices");
    }
    require(mesh.source.facePropertyIds.size() == mesh.faces.size() &&
        mesh.source.faceSourceLines.size() == mesh.faces.size(), "source table alignment");
}

void generate(Mesh& mesh, std::size_t components, std::size_t count, bool closed = true) {
    RWGOptions options;
    options.requireClosedSurface = closed;
    RWGGenerationReport report;
    require(RWGGenerator::GenerateRWG(mesh.faces, mesh.points, mesh.bases, options, report), "generate");
    require(report.success && report.components.size() == components && mesh.bases.size() == count,
        "analytic component/RWG count");
    require(report.rwgComponentIds.size() == count, "RWG component table alignment");
    std::array<int, 2> previous{0, 0};
    std::set<const FaceElement*> owned;
    for (const auto& face : mesh.faces) owned.insert(&face);
    for (std::size_t i = 0; i < mesh.bases.size(); ++i) {
        const auto& b = mesh.bases[i];
        const std::array<int, 2> key{b.commonEdgeVertex1->id, b.commonEdgeVertex2->id};
        require(key > previous && key[0] < key[1] && b.rwgid == static_cast<int>(i),
            "deterministic original-node edge ordering and global IDs");
        previous = key;
        require(owned.count(b.positiveFace) && owned.count(b.negativeFace), "RWG owns no temporary faces");
        require(b.edgeLength > 0 && std::isfinite(b.edgeLength), "valid edge length");
        require(report.faceComponentIds[static_cast<std::size_t>(b.positiveFace - mesh.faces.data())] ==
            report.faceComponentIds[static_cast<std::size_t>(b.negativeFace - mesh.faces.data())],
            "both support faces have the same component");
    }
}

static bool cubeNear(const OCTree::Node& a, const OCTree::Node& b) {
    const double limit = 1.01 * a.cubeLength;
    return std::fabs(a.center.x - b.center.x) <= limit &&
        std::fabs(a.center.y - b.center.y) <= limit &&
        std::fabs(a.center.z - b.center.z) <= limit;
}

static void verifyOctree() {
    const std::array<double, 4> centers{0.0, 1.0, 3.0, 1.0};
    std::vector<RWGBase> bases(centers.size());
    for (std::size_t i = 0; i < centers.size(); ++i) {
        bases[i].rwgid = static_cast<int>(i);
        bases[i].edgeCenterCenter = {0, centers[i], 0.0, 0.0};
    }

    const OCTree tree(bases, 2.0, "octree_fixture.nas");
    bool hasNear = false, hasDistant = false;
    for (int level = 2; level <= tree.maxLevel_; ++level) {
        const auto& nodes = tree.octreeNodes_[level];
        const auto& vectors = tree.octreeNodesDRvec_[level - 2];
        for (const auto* a : nodes) {
            std::set<const OCTree::Node*> near, distant;
            for (const auto* b : nodes) {
                if (a == b) continue;
                if (cubeNear(*a, *b)) near.insert(b);
                else if (a->realFather && b->realFather && cubeNear(*a->realFather, *b->realFather)) {
                    distant.insert(b);
                }
            }
            require(near == std::set<const OCTree::Node*>(a->nearNeighbors.begin(), a->nearNeighbors.end()),
                "octree near neighbors");
            require(distant == std::set<const OCTree::Node*>(a->distantRelatives.begin(), a->distantRelatives.end()),
                "octree distant relatives");
            hasNear = hasNear || !near.empty();
            hasDistant = hasDistant || !distant.empty();

            require(a->dRVecId.size() == a->distantRelatives.size(), "octree displacement count");
            for (std::size_t i = 0; i < a->dRVecId.size(); ++i) {
                const auto* b = a->distantRelatives[i];
                const auto id = static_cast<std::size_t>(a->dRVecId[i]);
                require(id < vectors.size(), "octree displacement index");
                const auto& offset = vectors[id];
                require(std::fabs(offset.x - (a->center.x - b->center.x)) <= 1e-12 &&
                    std::fabs(offset.y - (a->center.y - b->center.y)) <= 1e-12 &&
                    std::fabs(offset.z - (a->center.z - b->center.z)) <= 1e-12,
                    "octree displacement matches its distant relative");
            }
        }
    }
    require(hasNear && hasDistant, "octree fixture exercises near and distant interactions");

    std::vector<int> coverage(bases.size(), 0);
    for (const auto* leaf : tree.octreeNodes_[tree.maxLevel_]) {
        for (const auto* basis : leaf->rwgIndices) {
            const auto id = static_cast<std::size_t>(basis->rwgid);
            require(id < bases.size() && basis == &bases[id], "octree leaf RWG identity");
            ++coverage[id];
        }
    }
    require(std::all_of(coverage.begin(), coverage.end(), [](int count) { return count == 1; }),
        "every RWG belongs to exactly one leaf");
}

int main(int argc, char** argv) {
    require(argc == 3, "usage: nas_rwg_validation fixture_dir data_root");
    verifyOctree();
    const fs::path dir(argv[1]);
    const fs::path dataRoot(argv[2]);
    fs::create_directories(dir);
    const auto p1 = tetraPoints(0, 0);
    const auto f1 = tetraFaces(0, 1);
    const auto p2 = tetraPoints(100, 3);
    const auto f2 = tetraFaces(100, 11);
    const std::array<std::string, 3> arrangements{p1 + p2 + f1 + f2,
        p1 + f1 + p2 + f2, f1 + f2 + p1 + p2};
    std::vector<std::array<int, 4>> reference;
    for (std::size_t order = 0; order < arrangements.size(); ++order) {
        Mesh mesh;
        readMesh(mesh, writeFixture(dir, "ordering_" + std::to_string(order), arrangements[order]));
        generate(mesh, 2, 12);
        require(mesh.source.facePropertyIds[0] == 100 && mesh.source.facePropertyIds[3] == 103,
            "one component may span multiple PIDs");
        std::vector<std::array<int, 4>> signature;
        for (const auto& b : mesh.bases) signature.push_back({b.commonEdgeVertex1->id,
            b.commonEdgeVertex2->id, b.positiveFace->faceId, b.negativeFace->faceId});
        if (order == 0) reference = signature;
        else require(signature == reference, "record order does not alter supports or signs");
    }

    std::string many;
    for (int i = 0; i < 64; ++i) many += tetraPoints(i * 10, i * 3.0) + tetraFaces(i * 10, i * 4 + 1);
    Mesh growth;
    readMesh(growth, writeFixture(dir, "many_interleaved", many));
    generate(growth, 64, 384);

    Mesh large;
    const auto largeCards = largePoint(1, 0, 0, 0, "A") + largePoint(2, 1, 0, 0, "B") +
        largePoint(3, 0, 1, 0, "C") + largePoint(4, 0, 0, 1, "D");
    readMesh(large, writeFixture(dir, "large_fields", largeCards + f1));
    generate(large, 1, 6);

    Mesh small;
    std::ostringstream smallCards;
    const double xyz[4][3] = {{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
    for (int i = 0; i < 4; ++i) {
        smallCards << std::left << std::setw(8) << "GRID" << std::right
            << std::setw(8) << i + 1 << std::setw(8) << "";
        for (double v : xyz[i]) smallCards << std::setw(8) << v;
        smallCards << '\n';
    }
    readMesh(small, writeFixture(dir, "small_fields", smallCards.str() + f1));
    generate(small, 1, 6);

    Mesh exponents;
    readMesh(exponents, writeFixture(dir, "exponents",
        "GRID,1,,0,0,0\nGRID,2,,1.0D+0,0,0\nGRID,3,,0,1.0+0,0\nGRID,4,,0,0,+1.0E+0\n" + f1));
    generate(exponents, 1, 6);

    Mesh stable;
    readMesh(stable, writeFixture(dir, "stable", p1 + f1));
    generate(stable, 1, 6);
    Point* const oldPoints = stable.points.data();
    FaceElement* const oldFaces = stable.faces.data();
    const auto oldPids = stable.source.facePropertyIds;
    const std::vector<std::string> invalid{
        p1 + "GRID,1,,4,0,0\n" + f1,
        p1 + f1 + "CTRIA3,1,1,1,2,3\n",
        p1 + "CTRIA3,1,1,1,2,99\n",
        p1 + "CTRIA3,1,1,1,1,3\n",
        "GRID,1,,0,0,0\nGRID,2,,1,0,0\nGRID,3,,2,0,0\nCTRIA3,1,1,1,2,3\n",
        "GRID,1,,NaN,0,0\n" + f1,
        "GRID,1,5,0,0,0\n" + f1,
        "GRID,2147483648,,0,0,0\n" + f1,
        largeCards.substr(0, largeCards.find('\n') + 1),
        "CQUAD4,1,1,1,2,3,4\n",
        "GRID,1,,0.0junk,0,0\n" + f1,
        "$ empty mesh\n"};
    for (std::size_t i = 0; i < invalid.size(); ++i) {
        require(readNasData(writeFixture(dir, "invalid_" + std::to_string(i), invalid[i]).string(),
            stable.points, stable.faces, stable.source) != 0, "reject invalid input");
        require(stable.points.data() == oldPoints && stable.faces.data() == oldFaces &&
            stable.source.facePropertyIds == oldPids && stable.bases[0].positiveFace->triangularArea > 0,
            "failed read preserves original storage, metadata, and borrowed RWGs");
    }

    RWGOptions closed;
    closed.requireClosedSurface = true;
    RWGGenerationReport report;
    RWGBase* const oldBases = stable.bases.data();
    Mesh nonmanifold;
    readMesh(nonmanifold, writeFixture(dir, "nonmanifold", p1 + "GRID,5,,0,-1,0\n"
        "CTRIA3,20,1,1,2,3\nCTRIA3,21,1,2,1,4\nCTRIA3,22,1,1,2,5\n"));
    require(!RWGGenerator::GenerateRWG(nonmanifold.faces, nonmanifold.points,
        stable.bases, closed, report), "reject triple edge");
    require(!report.success && report.nonManifoldEdgeCount == 1 && report.components.size() == 1 &&
        report.nonManifoldEdges[0].faceIndices.size() == 3 && report.rwgComponentIds.empty(),
        "failure diagnostics retain every incident face");
    require(stable.bases.data() == oldBases && stable.bases.size() == 6, "failed generation preserves RWGs");

    Mesh open;
    readMesh(open, writeFixture(dir, "open", "GRID,1,,0,0,0\nGRID,2,,1,0,0\n"
        "GRID,3,,1,1,0\nGRID,4,,0,1,0\nCTRIA3,1,1,1,2,3\nCTRIA3,2,1,1,3,4\n"));
    generate(open, 1, 1, false);
    require(!RWGGenerator::GenerateRWG(open.faces, open.points, open.bases, closed, report) &&
        report.boundaryEdgeCount == 4 && open.bases.size() == 1, "closed policy rejects open mesh");
    require(RWGGenerator::GenerateRWG(open.faces, open.points, open.bases), "legacy default allows boundaries");

    Mesh duplicate;
    readMesh(duplicate, writeFixture(dir, "duplicate_triangle", p1 +
        "CTRIA3,1,1,1,2,3\nCTRIA3,2,1,3,2,1\n"));
    require(!RWGGenerator::GenerateRWG(duplicate.faces, duplicate.points, duplicate.bases), "duplicate geometry");
    Mesh isolated;
    readMesh(isolated, writeFixture(dir, "isolated", p1 + "CTRIA3,1,1,1,2,3\n"));
    require(!RWGGenerator::GenerateRWG(isolated.faces, isolated.points, isolated.bases), "zero RWG rejected");

    Mesh badWinding;
    std::string reversed = f1;
    reversed.replace(0, reversed.find('\n'), "CTRIA3,1,100,1,2,3");
    readMesh(badWinding, writeFixture(dir, "winding", p1 + reversed));
    require(!RWGGenerator::GenerateRWG(badWinding.faces, badWinding.points,
        badWinding.bases, closed, report) && report.orientationConflictCount == 3, "closed winding validation");

    Mesh pinched;
    const std::string pinchedPoints = p1 + "GRID,5,,-1,0,0\nGRID,6,,0,-1,0\nGRID,7,,0,0,-1\n";
    readMesh(pinched, writeFixture(dir, "pinched_vertex", pinchedPoints + f1 +
        "CTRIA3,11,1,1,6,5\nCTRIA3,12,1,1,5,7\nCTRIA3,13,1,5,6,7\nCTRIA3,14,1,6,1,7\n"));
    require(!RWGGenerator::GenerateRWG(pinched.faces, pinched.points, pinched.bases, closed, report) &&
        report.nonManifoldEdgeCount == 0 && report.error.find("fans") != std::string::npos,
        "closed vertex-link validation");

    Mesh cube;
    readMesh(cube, dataRoot / "Cube_0d5m_Die_1e9/Cube_0d5m_Die_1e9.nas");
    generate(cube, 1, 3402);
    double area = 0;
    for (const auto& face : cube.faces) area += face.triangularArea;
    require(std::abs(area - 1.5) < 1e-8, "cube total area");
    Mesh sphere;
    readMesh(sphere, dataRoot / "Sphere_1e9/Sphere_1e9.nas");
    generate(sphere, 1, 8127);
    Mesh four;
    readMesh(four, dataRoot / "4Qiuzhui/4Qiuzhui_2e9.nas");
    require(four.points.size() == 4875 && four.faces.size() == 9886, "four-target input fully read");
    require(!RWGGenerator::GenerateRWG(four.faces, four.points, four.bases, closed, report) &&
        report.components.size() == 4 && report.nonManifoldEdgeCount == 160 &&
        report.boundaryEdgeCount == 0, "four-target topology independently audited counts");
    require(four.source.faceSourceLines[0] == 2452 && report.nonManifoldEdges[0].vertexIds ==
        std::array<int, 2>{2, 3}, "real-data diagnostic source mapping");

    std::cout << "PASS: NAS/RWG ordering, pointer lifetime, rollback, topology, formats, and real meshes\n";
    return 0;
}
