// NasReader.cpp
#include "NASReader.h"
#include "PathUtils.h"

#include <fstream>
#include <iostream>
#include <iomanip>

using namespace std;

Point FindPoint(int id, const std::vector<Point>& points) {
	for (const auto& p : points) {
		if (p.id == id) return p;
	}
	throw std::runtime_error("Point ID " + std::to_string(id) + " not found");
}

Point CalculateTriangleCenter(const Point& a, const Point& b, const Point& c) {
	return { 0,
			(a.x + b.x + c.x) / 3.0,
			(a.y + b.y + c.y) / 3.0,
			(a.z + b.z + c.z) / 3.0 };
}

static double CalculateTriangleArea(const Point& a, const Point& b, const Point& c) {
	Point ab = b - a;
	Point ac = c - a;

	Point cross = {
		0,
		ab.y * ac.z - ab.z * ac.y,
		ab.z * ac.x - ab.x * ac.z,
		ab.x * ac.y - ab.y * ac.x
	};

	return sqrt(cross.x * cross.x + cross.y * cross.y + cross.z * cross.z) / 2.0;
}

static Point CalculateNormalVector(const Point& a, const Point& b, const Point& c) {
	Point ab = b - a;
	Point ac = c - a;

	Point normal = {
		0,
		ab.y * ac.z - ab.z * ac.y,
		ab.z * ac.x - ab.x * ac.z,
		ab.x * ac.y - ab.y * ac.x
	};

	double length = sqrt(normal.x * normal.x + normal.y * normal.y + normal.z * normal.z);
	if (length > 0) {
		normal = normal / length;
	}
	if (std::abs(normal.x) < 1e-6)  normal.x = 0;
	if (std::abs(normal.y) < 1e-6)  normal.y = 0;
	if (std::abs(normal.z) < 1e-6)  normal.z = 0;
	return normal;
}

int readNasData(const string& filename, vector<Point>& points, vector<FaceElement>& faces) {
	ifstream file(filename);
	if (!file.is_open()) {
		cerr << "Error opening file: " << filename << endl;
		return -1;
	}

	string line;
	bool expectingGridContinuation = false;
	Point tempPoint;

	while (getline(file, line)) {
		if (line.substr(0, 5) == "GRID*") {
			istringstream iss(line.substr(5));
			if (!(iss >> tempPoint.id >> tempPoint.x >> tempPoint.y)) {
				cerr << "Invalid GRID* format: " << line << endl;
				return -1;
			}
			expectingGridContinuation = true;
		}
		else if (expectingGridContinuation && line[0] == '*') {
			istringstream iss(line.substr(1));
			int continuationId;
			double zValue;
			if (!(iss >> continuationId >> zValue)) {
				cerr << "Invalid continuation line: " << line << endl;
				return -1;
			}

			if (continuationId != tempPoint.id) {
				cerr << "ID mismatch between GRID* and * lines" << endl;
				return -1;
			}

			tempPoint.z = zValue;
			points.push_back(tempPoint);
			expectingGridContinuation = false;
		}
		else if (line.substr(0, 6) == "CTRIA3") {
			FaceElement element;
			int propertyId;
			istringstream iss(line.substr(6));

			int vertex1Id, vertex2Id, vertex3Id;
			if (!(iss >> element.faceId >> propertyId
				>> vertex1Id >> vertex2Id >> vertex3Id)) {
				cerr << "Invalid CTRIA3 format: " << line << endl;
				return -1;
			}

			for (auto& p : points) {
				if (p.id == vertex1Id) element.vertex1 = &p;
				if (p.id == vertex2Id) element.vertex2 = &p;
				if (p.id == vertex3Id) element.vertex3 = &p;
			}
			if (!element.vertex1 || !element.vertex2 || !element.vertex3) {
				std::cerr << "Error: Vertex not found for face " << element.faceId << " in " << filename << std::endl;
				return -1;
			}
			faces.push_back(element);
		}
	}

	for (auto& face : faces) {
		if (!face.vertex1 || !face.vertex2 || !face.vertex3) {
			std::cerr << "Error: Null vertex pointer in face " << face.faceId << std::endl;
			return -1;
		}

		Point p1 = *face.vertex1;
		Point p2 = *face.vertex2;
		Point p3 = *face.vertex3;

		face.triangularArea = CalculateTriangleArea(p1, p2, p3);
		face.triangleCenter = CalculateTriangleCenter(p1, p2, p3);
		face.externalNormalVector = CalculateNormalVector(p1, p2, p3);
	}

	if (expectingGridContinuation) {
		cerr << "Unfinished GRID* entry at end of file" << endl;
		return -1;
	}
	cout << "Info: NAS file read (NASReader.cpp)" << endl;
	return 0;
}

bool writeOutputFiles(const string& nasPath, const vector<FaceElement>& faces, const vector<Point>& points) {
	string parentPath = getParentPath(nasPath);
	string baseName = getBaseName(nasPath);

	ofstream faceFile(parentPath + baseName + "_DEBUG_Faces.txt");
	if (!faceFile) {
		cerr << "Error creating faces output file" << endl;
		return false;
	}

	faceFile << "FaceID\tVertex1\tVertex2\tVertex3\tArea\tCenterX\tCenterY\tCenterZ\tExternalNormalVector(x)\tExternalNormalVector(y)\tExternalNormalVector(z)\n";
	for (const auto& face : faces) {
		faceFile << face.faceId << "\t"
			<< face.vertex1->id << "\t"
			<< face.vertex2->id << "\t"
			<< face.vertex3->id << "\t"
			<< face.triangularArea << "\t"
			<< face.triangleCenter.x << "\t"
			<< face.triangleCenter.y << "\t"
			<< face.triangleCenter.z << "\t"
			<< face.externalNormalVector.x << "\t"
			<< face.externalNormalVector.y << "\t"
			<< face.externalNormalVector.z << "\n";
	}

	cout << "Info: Face element data written to: " << parentPath + baseName + "_DEBUG_Faces.txt (NASReader.cpp)" << std::endl;

	ofstream pointFile(parentPath + baseName + "_DEBUG_Points.txt");
	if (!pointFile) {
		cerr << "Error creating points output file" << endl;
		return false;
	}

	pointFile << fixed << setprecision(6);
	pointFile << "PointID\tX\tY\tZ\n";
	for (const auto& point : points) {
		pointFile << point.id << "\t"
			<< point.x << "\t"
			<< point.y << "\t"
			<< point.z << "\n";
	}
	cout << "Info: Point data written to: " << parentPath + baseName + "_DEBUG_Points.txt (NASReader.cpp)" << std::endl;

	return true;
}