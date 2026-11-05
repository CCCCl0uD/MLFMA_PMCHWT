// RWGGenerator.cpp
#include "RWGGenerator.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <numeric>
#include <set>
#include <unordered_map>
#include <unordered_set>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <sstream>

using namespace std;

bool RWGGenerator::GenerateRWG(const vector<FaceElement>& faces,
    const vector<Point>& points, vector<RWGBase>& result) {
    RWGGenerationReport report;
    return GenerateRWG(faces, points, result, RWGOptions{}, report);
}

bool RWGGenerator::GenerateRWG(const vector<FaceElement>& faces,
    const vector<Point>&, vector<RWGBase>& result,
    const RWGOptions& options, RWGGenerationReport& report) {
    RWGGenerationReport analysis;
    const auto fail = [&](string message) {
        analysis.error = move(message);
        analysis.rwgComponentIds.clear();
        cerr << "Error: " << analysis.error << " (RWGGenerator.cpp)\n";
        report.swap(analysis);
        return false;
    };
    // Input geometry was validated by readNasData. Only topology is checked here.
    struct EdgeUse {
        size_t faceIndex;
        size_t localEdge;
        Point* freeVertex;
        bool ascending;
    };
    struct EdgeData {
        const Point* lowVertex;
        const Point* highVertex;
        vector<EdgeUse> uses;
    };
    map<Edge, EdgeData> edgeMap;
    set<array<int, 3>> triangles;
    unordered_map<int, vector<size_t>> vertexFaces;
    for (size_t i = 0; i < faces.size(); ++i) {
        const auto& face = faces[i];
        Point* vertices[3] = {face.vertex1, face.vertex2, face.vertex3};
        array<int, 3> triangle{vertices[0]->id, vertices[1]->id, vertices[2]->id};
        sort(triangle.begin(), triangle.end());
        if (!triangles.insert(triangle).second)
            return fail("Duplicate triangle at EID " + to_string(face.faceId));
        for (size_t j = 0; j < 3; ++j) {
            Point* a = vertices[j];
            Point* b = vertices[(j + 1) % 3];
            const bool ascending = a->id < b->id;
            auto inserted = edgeMap.try_emplace(Edge(a->id, b->id),
                EdgeData{ascending ? a : b, ascending ? b : a, {}});
            inserted.first->second.uses.push_back({i, j, vertices[(j + 2) % 3], ascending});
            vertexFaces[a->id].push_back(i);
        }
    }

    vector<size_t> parents(faces.size());
    iota(parents.begin(), parents.end(), size_t{0});
    const auto findRoot = [&](size_t face) {
        while (parents[face] != face) {
            parents[face] = parents[parents[face]];
            face = parents[face];
        }
        return face;
    };
    vector<array<size_t, 3>> neighbors(faces.size(),
        array<size_t, 3>{faces.size(), faces.size(), faces.size()});
    for (const auto& entry : edgeMap) {
        const auto& uses = entry.second.uses;
        const size_t root = findRoot(uses[0].faceIndex);
        for (size_t j = 1; j < uses.size(); ++j) parents[findRoot(uses[j].faceIndex)] = root;
        if (uses.size() == 1) ++analysis.boundaryEdgeCount;
        else if (uses.size() == 2) {
            neighbors[uses[0].faceIndex][uses[0].localEdge] = uses[1].faceIndex;
            neighbors[uses[1].faceIndex][uses[1].localEdge] = uses[0].faceIndex;
            if (uses[0].ascending == uses[1].ascending) ++analysis.orientationConflictCount;
        }
        else {
            ++analysis.nonManifoldEdgeCount;
            RWGEdgeIssue issue;
            issue.vertexIds = {entry.first.v[0], entry.first.v[1]};
            for (const auto& use : uses) issue.faceIndices.push_back(use.faceIndex);
            analysis.nonManifoldEdges.push_back(move(issue));
        }
    }

    unordered_map<size_t, size_t> componentIndices;
    analysis.faceComponentIds.resize(faces.size());
    for (size_t i = 0; i < faces.size(); ++i) {
        const auto inserted = componentIndices.emplace(findRoot(i), componentIndices.size());
        if (inserted.second) analysis.components.emplace_back();
        const size_t component = inserted.first->second;
        analysis.faceComponentIds[i] = component;
        ++analysis.components[component].faceCount;
    }
    vector<unordered_set<int>> componentVertices(analysis.components.size());
    for (size_t i = 0; i < faces.size(); ++i) {
        auto& ids = componentVertices[analysis.faceComponentIds[i]];
        ids.insert(faces[i].vertex1->id);
        ids.insert(faces[i].vertex2->id);
        ids.insert(faces[i].vertex3->id);
    }
    for (size_t i = 0; i < analysis.components.size(); ++i)
        analysis.components[i].vertexCount = componentVertices[i].size();
    size_t rwgCount = 0;
    for (const auto& entry : edgeMap) {
        const auto& uses = entry.second.uses;
        auto& component = analysis.components[analysis.faceComponentIds[uses[0].faceIndex]];
        ++component.edgeCount;
        if (uses.size() == 2) {
            ++component.rwgCount;
            ++rwgCount;
        }
    }

    cout << "Info: Mesh components: " << analysis.components.size()
        << ", boundary edges: " << analysis.boundaryEdgeCount
        << ", non-manifold edges: " << analysis.nonManifoldEdgeCount
        << ", orientation conflicts: " << analysis.orientationConflictCount << '\n';
    if (analysis.nonManifoldEdgeCount != 0) {
        for (const auto& issue : analysis.nonManifoldEdges) {
            cerr << "  Edge (" << issue.vertexIds[0] << ", " << issue.vertexIds[1] << "), EIDs:";
            for (size_t i : issue.faceIndices) cerr << ' ' << faces[i].faceId;
            cerr << '\n';
        }
        return fail("Non-manifold edges are not supported by ordinary RWGs");
    }
    if (options.requireClosedSurface && analysis.boundaryEdgeCount != 0)
        return fail("This integral equation requires a closed surface");
    if (options.requireClosedSurface && analysis.orientationConflictCount != 0)
        return fail("Inconsistent triangle winding on a closed surface");
    if (rwgCount == 0) return fail("No interior edges: no usable RWG unknowns");

    if (options.requireClosedSurface) {
        // A closed vertex link must be one connected fan, not several shells
        // sharing only a node. Edge multiplicity was checked above.
        vector<size_t> visited(faces.size(), 0);
        size_t vertexStamp = 0;
        for (const auto& vertex : vertexFaces) {
            ++vertexStamp;
            const auto& incident = vertex.second;
            vector<size_t> pending{incident[0]};
            visited[incident[0]] = vertexStamp;
            size_t reached = 0;
            while (!pending.empty()) {
                const size_t i = pending.back();
                pending.pop_back();
                ++reached;
                const auto& face = faces[i];
                Point* vertices[3] = {face.vertex1, face.vertex2, face.vertex3};
                for (size_t j = 0; j < 3; ++j) {
                    if (vertices[j]->id != vertex.first && vertices[(j + 1) % 3]->id != vertex.first) continue;
                    const size_t next = neighbors[i][j];
                    if (visited[next] != vertexStamp) {
                        visited[next] = vertexStamp;
                        pending.push_back(next);
                    }
                }
            }
            if (reached != incident.size())
                return fail("Disconnected triangle fans at GRID " + to_string(vertex.first));
        }
    }

    vector<RWGBase> bases;
    bases.reserve(rwgCount);
    analysis.rwgComponentIds.reserve(rwgCount);
    for (const auto& entry : edgeMap) {
        const auto& edge = entry.second;
        if (edge.uses.size() == 1) continue; // Explicitly permitted boundary edges.
        const auto& positive = edge.uses[0];
        const auto& negative = edge.uses[1];
        RWGBase base{};
        base.rwgid = static_cast<int>(bases.size());
        base.positiveFace = const_cast<FaceElement*>(&faces[positive.faceIndex]);
        base.negativeFace = const_cast<FaceElement*>(&faces[negative.faceIndex]);
        base.commonEdgeVertex1 = edge.lowVertex;
        base.commonEdgeVertex2 = edge.highVertex;
        base.freeVertexPositive = positive.freeVertex;
        base.freeVertexNegative = negative.freeVertex;
        const Point delta = *edge.highVertex - *edge.lowVertex;
        base.edgeLength = hypot(delta.x, delta.y, delta.z);
        base.edgeCenterCenter = {0,
            (edge.lowVertex->x + edge.highVertex->x) / 2.0,
            (edge.lowVertex->y + edge.highVertex->y) / 2.0,
            (edge.lowVertex->z + edge.highVertex->z) / 2.0};
        if (!(base.edgeLength > 0.0) || !isfinite(base.edgeLength) ||
            !isfinite(base.edgeCenterCenter.x) || !isfinite(base.edgeCenterCenter.y) ||
            !isfinite(base.edgeCenterCenter.z)) return fail("Invalid common-edge geometry");
        bases.push_back(base);
        analysis.rwgComponentIds.push_back(analysis.faceComponentIds[positive.faceIndex]);
    }
    analysis.success = true;
    cout << "Info: RWG basis generation done, total: " << bases.size()
        << " (RWGGenerator.cpp)\n";
    result.swap(bases);
    report.swap(analysis);
    return true;
}

bool RWGGenerator::WriteRWGToFile(const vector<RWGBase>& bases, const string& nasPath) {
    const auto separator = nasPath.find_last_of("/\\");
    const auto filenameStart = separator == string::npos ? 0 : separator + 1;
    const auto extension = nasPath.find_last_of('.');
    const string outputPath = (extension != string::npos && extension > filenameStart
        ? nasPath.substr(0, extension) : nasPath) + "_DEBUG_RWG.txt";
    ofstream outFile(outputPath);
    if (!outFile) {
        cerr << "Error: Cannot create output file: " << outputPath << '\n';
        return false;
    }
    outFile << "RWGID\tPositiveFace\tNegativeFace\t"
        << "CommonEdgeVertex1\tCommonEdgeVertex2\t"
        << "FreeVertexPos\tFreeVertexNeg\tEdgeLength\t"
        << "EdgeCenter_X\tEdgeCenter_Y\tEdgeCenter_Z\n";
    outFile << fixed << setprecision(6);
    for (const auto& b : bases) {
        outFile << b.rwgid << '\t' << b.positiveFace->faceId << '\t'
            << b.negativeFace->faceId << '\t' << b.commonEdgeVertex1->id << '\t'
            << b.commonEdgeVertex2->id << '\t' << b.freeVertexPositive->id << '\t'
            << b.freeVertexNegative->id << '\t' << b.edgeLength << '\t'
            << b.edgeCenterCenter.x << '\t' << b.edgeCenterCenter.y << '\t'
            << b.edgeCenterCenter.z << '\n';
    }
    outFile.close();
    if (!outFile) {
        cerr << "Error writing RWG output file: " << outputPath << '\n';
        return false;
    }
    cout << "Info: RWG basis written to: " << outputPath << " (RWGGenerator.cpp)\n";
    return true;
}
