// MyStruct.h
#pragma once
#ifndef MYSTRUCT_H
#define MYSTRUCT_H

#include <array>
#include <complex>

// ====== Complex3D ======
struct Complex3D {
	std::complex<double> x;
	std::complex<double> y;
	std::complex<double> z;

	Complex3D();
	Complex3D(const std::complex<double>& x_, const std::complex<double>& y_, const std::complex<double>& z_);

	Complex3D operator+(const Complex3D& other) const;
	Complex3D operator*(double scalar) const;
	Complex3D operator*(const std::complex<double>& scalar) const;
};

// ====== kp_Point ======
struct kp_Point {
	double x;
	double y;
	double z;

	kp_Point() : x(0), y(0), z(0) {}
	kp_Point(double x_, double y_, double z_) : x(x_), y(y_), z(z_) {}
};

// ====== nodePoint ======
struct nodePoint {
	double x, y, z;
	nodePoint(double x = 0, double y = 0, double z = 0) : x(x), y(y), z(z) {}
};

// ====== InterpData ======
struct InterpData {
	int fIdx = -1;
	std::array<int, 3> tIdx{ -1, -1, -1 };
	std::array<int, 3> pIdx{ -1, -1, -1 };
	std::array<int, 9> tp{ -1, -1, -1, -1, -1, -1, -1, -1, -1 };
	std::array<double, 9> w{};

	InterpData() = default;
	InterpData(int f, const std::array<int, 3>& t, const std::array<int, 3>& p,
		const std::array<double, 9>& weight, int phiCount)
		: fIdx(f), tIdx(t), pIdx(p), w(weight) {
		int k = 0;
		for (int i = 0; i < 3; ++i)
			for (int j = 0; j < 3; ++j)
				tp[k++] = tIdx[i] * phiCount + pIdx[j];
	}
};

// ====== Point ======
struct Point {
	int id;
	double x;
	double y;
	double z;

	Point operator+(const Point& p) const {
		return { 0, x + p.x, y + p.y, z + p.z };
	}
	Point operator-(const Point& p) const {
		return { 0, x - p.x, y - p.y, z - p.z };
	}
	Point operator*(double scalar) const {
		return { 0, x * scalar, y * scalar, z * scalar };
	}
	Point operator/(double scalar) const {
		return { 0, x / scalar, y / scalar, z / scalar };
	}
};

// ====== FaceElement ======
struct FaceElement {
	int faceId;
	Point* vertex1 = nullptr;
	Point* vertex2 = nullptr;
	Point* vertex3 = nullptr;
	Point externalNormalVector{};
	double triangularArea = 0.0;
	Point triangleCenter{};
};

// ====== RWGBase ======
struct RWGBase {
	int rwgid;
	FaceElement* positiveFace;
	FaceElement* negativeFace;
	const Point* commonEdgeVertex1;
	const Point* commonEdgeVertex2;
	Point* freeVertexPositive;
	Point* freeVertexNegative;
	double edgeLength;
	Point edgeCenterCenter;
};

#endif // MYSTRUCT_H