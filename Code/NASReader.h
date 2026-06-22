// NasReader.h
#pragma once
#ifndef NASREADER_H
#define NASREADER_H

#include "MyStruct.h"

#include <string>
#include <vector>

Point FindPoint(int id, const std::vector<Point>& points);
Point CalculateTriangleCenter(const Point& a, const Point& b, const Point& c);
int readNasData(const std::string& filename, std::vector<Point>& points, std::vector<FaceElement>& faces);
bool writeOutputFiles(const std::string& nasPath, const std::vector<FaceElement>& faces, const std::vector<Point>& points);

#endif //NASREADER_H