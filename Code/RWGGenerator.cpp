// RWGGenerator.cpp
#include "RWGGenerator.h"

#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <sstream>

using namespace std;

namespace PathUtil {
	string GetParentPath(const string& fullPath) {
		size_t pos = fullPath.find_last_of("/\\");
		return (pos != string::npos) ? fullPath.substr(0, pos + 1) : "";
	}

	string GetBaseName(const string& fullPath) {
		size_t start = fullPath.find_last_of("/\\");
		start = (start == string::npos) ? 0 : start + 1;
		size_t end = fullPath.find_last_of('.');
		return (end != string::npos && end > start) ?
			fullPath.substr(start, end - start) : fullPath.substr(start);
	}
}

static Point FindPoint(int id, const vector<Point>& points) {
	for (const auto& p : points) {
		if (p.id == id) return p;
	}
	throw runtime_error("Point ID " + to_string(id) + " not found");
}

static double CalcDistance(const Point& a, const Point& b) {
	return sqrt(pow(a.x - b.x, 2) + pow(a.y - b.y, 2) + pow(a.z - b.z, 2));
}

bool RWGGenerator::GenerateRWG(const vector<FaceElement>& faces, const vector<Point>& points, vector<RWGBase>& result) {
	try {
		map<Edge, vector<size_t>> edgeMap;

		for (size_t i = 0; i < faces.size(); ++i) {
			const auto& f = faces[i];

			Edge edges[3] = {
				Edge(f.vertex1->id, f.vertex2->id),
				Edge(f.vertex2->id, f.vertex3->id),
				Edge(f.vertex3->id, f.vertex1->id)
			};

			for (int j = 0; j < 3; ++j) {
				edgeMap[edges[j]].push_back(i);
			}
		}

		vector<RWGBase> bases;
		size_t rwg_id_counter = 0;

		for (const auto& entry : edgeMap) {
			const auto& edge = entry.first;
			const auto& faceIndices = entry.second;

			if (faceIndices.size() != 2) continue;

			const size_t posIdx = faceIndices[0];
			const size_t negIdx = faceIndices[1];
			const auto& facePos = &faces[posIdx];
			const auto& faceNeg = &faces[negIdx];

			RWGBase base;
			base.rwgid = rwg_id_counter++;  // ����ΨһID
			base.positiveFace = const_cast<FaceElement*>(facePos);
			base.negativeFace = const_cast<FaceElement*>(faceNeg);

			for (auto& p : points) {
				if (p.id == edge.v[0]) base.commonEdgeVertex1 = &p;
				if (p.id == edge.v[1]) base.commonEdgeVertex2 = &p;
			}

			auto findFreeVertex = [&](const FaceElement& face) -> Point* {
				for (Point* vid : { face.vertex1, face.vertex2, face.vertex3 }) {
					if (vid->id != edge.v[0] && vid->id != edge.v[1]) {
						return vid;
					}
				}
				throw logic_error("Free vertex not found");
				};
			base.freeVertexPositive = findFreeVertex(*facePos);
			base.freeVertexNegative = findFreeVertex(*faceNeg);

			Point p1 = *base.commonEdgeVertex1;
			Point p2 = *base.commonEdgeVertex2;
			base.edgeLength = CalcDistance(p1, p2);
			base.edgeCenterCenter = {
				0,
				(p1.x + p2.x) / 2.0,
				(p1.y + p2.y) / 2.0,
				(p1.z + p2.z) / 2.0
			};

			bases.push_back(base);
		}

		result = move(bases);
		cout << "Info: RWG basis generation done, total: " << rwg_id_counter << " (RWGGenerator.cpp)" << endl;
		return true;
	}
	catch (const exception& e) {
		cerr << "[Error] " << e.what() << endl;
		return false;
	}
}

bool RWGGenerator::WriteRWGToFile(const vector<RWGBase>& bases, const string& nasPath) {
	try {
		const string parentPath = PathUtil::GetParentPath(nasPath);
		const string basePath = PathUtil::GetBaseName(nasPath);
		const string outputPath = parentPath + basePath + "_DEBUG_RWG.txt";

		ofstream outFile(outputPath);
		if (!outFile) {
			throw runtime_error("Cannot creat output file: " + outputPath);
		}

		outFile << "RWGID\tPositiveFace\tNegativeFace\t"
			<< "CommonEdgeVertex1\tCommonEdgeVertex2\t"
			<< "FreeVertexPos\tFreeVertexNeg\t"
			<< "EdgeLength\t"
			<< "EdgeCenter_X\tEdgeCenter_Y\tEdgeCenter_Z\n";

		outFile << fixed << setprecision(6);
		for (const auto& b : bases) {
			outFile << b.rwgid << "\t"
				<< b.positiveFace->faceId << "\t"
				<< b.negativeFace->faceId << "\t"
				<< b.commonEdgeVertex1->id << "\t"
				<< b.commonEdgeVertex2->id << "\t"
				<< b.freeVertexPositive->id << "\t"
				<< b.freeVertexNegative->id << "\t"
				<< b.edgeLength << "\t"
				<< b.edgeCenterCenter.x << "\t" << b.edgeCenterCenter.y << "\t" << b.edgeCenterCenter.z << "\n";
		}

		cout << "Info: RWG basis written to: " << outputPath << " (RWGGenerator.cpp)" << endl;
		return true;
	}
	catch (const exception& e) {
		cerr << "[Error] WriteRWGToFile: " << e.what() << endl;
		return false;
	}
}