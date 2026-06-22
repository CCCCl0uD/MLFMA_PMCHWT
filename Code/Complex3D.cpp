// Complex3D.cpp
#include "MyStruct.h"

Complex3D::Complex3D() : x(0.0, 0.0), y(0.0, 0.0), z(0.0, 0.0) {}
Complex3D::Complex3D(const std::complex<double>& x_, const std::complex<double>& y_, const std::complex<double>& z_)
	: x(x_), y(y_), z(z_) {
}

Complex3D Complex3D::operator+(const Complex3D& other) const {
	return Complex3D(x + other.x, y + other.y, z + other.z);
}
Complex3D Complex3D::operator*(double scalar) const {
	return Complex3D(x * scalar, y * scalar, z * scalar);
}
Complex3D Complex3D::operator*(const std::complex<double>& scalar) const {
	return Complex3D(x * scalar, y * scalar, z * scalar);
}