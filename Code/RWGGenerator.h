// RWGGenerator.h
#pragma once
#ifndef RWGGENERATOR_H
#define RWGGENERATOR_H
#include "NASReader.h"
#include <array>
#include <string>
#include <utility>
#include <vector>

struct RWGOptions {
    bool requireClosedSurface = false;
};

struct RWGComponentSummary {
    std::size_t vertexCount = 0;
    std::size_t faceCount = 0;
    std::size_t edgeCount = 0;
    std::size_t rwgCount = 0;
};

struct RWGEdgeIssue {
    std::array<int, 2> vertexIds{};
    // Index faces[]; source lines/PIDs are available in MeshReadReport.
    std::vector<std::size_t> faceIndices;
};

// Describes the latest attempt, including topology when generation fails.
// On failure success=false, rwgComponentIds is empty, and result is unchanged.
struct RWGGenerationReport {
    bool success = false;
    std::size_t boundaryEdgeCount = 0;
    std::size_t nonManifoldEdgeCount = 0;
    std::size_t orientationConflictCount = 0;
    std::vector<std::size_t> faceComponentIds;
    std::vector<std::size_t> rwgComponentIds;
    std::vector<RWGComponentSummary> components;
    std::vector<RWGEdgeIssue> nonManifoldEdges;
    std::string error;

    void swap(RWGGenerationReport& other) noexcept {
        using std::swap;
        swap(success, other.success);
        swap(boundaryEdgeCount, other.boundaryEdgeCount);
        swap(nonManifoldEdgeCount, other.nonManifoldEdgeCount);
        swap(orientationConflictCount, other.orientationConflictCount);
        faceComponentIds.swap(other.faceComponentIds);
        rwgComponentIds.swap(other.rwgComponentIds);
        components.swap(other.components);
        nonManifoldEdges.swap(other.nonManifoldEdges);
        error.swap(other.error);
    }
};

class RWGGenerator {
public:
    // Borrowed mesh: pass geometrically valid output from readNasData.
    // points/faces must stay alive and must not be resized, reordered or copied
    // while their RWGs are in use. Generation does not modify the mesh.
    static bool GenerateRWG(const std::vector<FaceElement>& faces,
        const std::vector<Point>& points, std::vector<RWGBase>& result);
    static bool GenerateRWG(const std::vector<FaceElement>& faces,
        const std::vector<Point>& points, std::vector<RWGBase>& result,
        const RWGOptions& options, RWGGenerationReport& report);

    static bool WriteRWGToFile(const std::vector<RWGBase>& bases,
        const std::string& nasPath);

private:
    struct Edge {
        int v[2];
        Edge(int a, int b) : v{(a < b) ? a : b, (a < b) ? b : a} {}
        bool operator<(const Edge& other) const {
            return (v[0] < other.v[0]) ||
                (v[0] == other.v[0] && v[1] < other.v[1]);
        }
    };
};

#endif // RWGGENERATOR_H
