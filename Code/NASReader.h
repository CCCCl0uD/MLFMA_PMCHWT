// NasReader.h
#pragma once
#ifndef NASREADER_H
#define NASREADER_H

#include "MyStruct.h"

#include <cstddef>
#include <string>
#include <vector>

// These arrays follow faces[] and belong to the same immutable mesh.
struct MeshReadReport {
	std::vector<int> facePropertyIds;
	std::vector<std::size_t> faceSourceLines;

	void swap(MeshReadReport& other) noexcept {
		facePropertyIds.swap(other.facePropertyIds);
		faceSourceLines.swap(other.faceSourceLines);
	}
};

Point CalculateTriangleCenter(const Point& a, const Point& b, const Point& c);
// Replaces outputs only on success. Do not replace a mesh borrowed by RWGs/tree.
int readNasData(const std::string& filename, std::vector<Point>& points, std::vector<FaceElement>& faces);
int readNasData(const std::string& filename, std::vector<Point>& points,
	std::vector<FaceElement>& faces, MeshReadReport& report);
bool writeOutputFiles(const std::string& nasPath, const std::vector<FaceElement>& faces, const std::vector<Point>& points);

#endif //NASREADER_H
