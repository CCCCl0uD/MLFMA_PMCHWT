// MovingTarget.h
#pragma once
#ifndef MOVING_TARGET_H
#define MOVING_TARGET_H

#include <array>
#include <string>
#include <vector>

struct MovingState {
	double time_s = 0.0;
	double position[3] = { 0.0, 0.0, 0.0 };
	double quat[4] = { 1.0, 0.0, 0.0, 0.0 };
};

struct MovingAngleSample {
	double time_s = 0.0;
	double pos[3] = { 0.0, 0.0, 0.0 };
	double range = 0.0;
	double thetaGlobalDeg = 0.0;
	double phiGlobalDeg = 0.0;
	double thetaBodyDeg = 0.0;
	double phiBodyDeg = 0.0;
};

namespace MovingTarget {
	std::vector<MovingState> readMovingPathFile(const std::string& path);
	void normalizeQuaternion(double q[4]);
	void quaternionToMatrix(const double q[4], double R[3][3]);
	std::array<double, 3> normalizeVector(const std::array<double, 3>& v);
	std::array<double, 3> rotateGlobalToBody(
		const double R[3][3],
		const std::array<double, 3>& vGlobal);
	void vectorToThetaPhiDeg(
		const std::array<double, 3>& v,
		double& thetaDeg,
		double& phiDeg);
	MovingAngleSample computeMovingAngleSample(
		const MovingState& state,
		const std::array<double, 3>& radarPos);
	std::vector<MovingAngleSample> computeMovingAngleSamples(
		const std::vector<MovingState>& states,
		const std::array<double, 3>& radarPos);
	std::vector<MovingAngleSample> mergeMovingAngleSamples(
		const std::vector<MovingAngleSample>& samples,
		double mergeTolDeg);
}

#endif // MOVING_TARGET_H
