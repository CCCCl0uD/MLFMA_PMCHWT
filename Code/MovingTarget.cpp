// MovingTarget.cpp
#include "MovingTarget.h"
#include "PhyConstant.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <sstream>
#include <stdexcept>

namespace {
	constexpr double kRangeEps = 1.0e-12;
	constexpr double kAnglePoleEpsDeg = 1.0e-9;

	std::string trimMovingLine(const std::string& s) {
		const size_t start = s.find_first_not_of(" \t\r\n");
		if (start == std::string::npos) {
			return "";
		}
		const size_t end = s.find_last_not_of(" \t\r\n");
		return s.substr(start, end - start + 1);
	}

	double normalizePhiDegLocal(double phiDeg) {
		double v = std::fmod(phiDeg, 360.0);
		if (v < 0.0) {
			v += 360.0;
		}
		return v;
	}

	bool isThetaPoleLocal(double thetaDeg) {
		return std::abs(thetaDeg) < kAnglePoleEpsDeg ||
			std::abs(thetaDeg - 180.0) < kAnglePoleEpsDeg;
	}

	double deltaPhiDeg(double a, double b) {
		const double da = normalizePhiDegLocal(a);
		const double db = normalizePhiDegLocal(b);
		const double d = std::abs(da - db);
		return std::min(d, 360.0 - d);
	}
}

namespace MovingTarget {
	std::vector<MovingState> readMovingPathFile(const std::string& path) {
		std::ifstream file(path);
		if (!file.is_open()) {
			throw std::runtime_error("Cannot open moving path file: " + path);
		}

		std::vector<MovingState> states;
		std::string line;
		int lineNo = 0;
		while (std::getline(file, line)) {
			++lineNo;
			const std::string trimmed = trimMovingLine(line);
			if (trimmed.empty() || trimmed[0] == '#') {
				continue;
			}

			std::istringstream iss(trimmed);
			MovingState state;
			if (!(iss >> state.time_s
				>> state.position[0] >> state.position[1] >> state.position[2]
				>> state.quat[0] >> state.quat[1] >> state.quat[2] >> state.quat[3])) {
				throw std::runtime_error(
					"Invalid moving path line " + std::to_string(lineNo) +
					": expected time x y z qw qx qy qz");
			}
			normalizeQuaternion(state.quat);
			states.push_back(state);
		}

		if (states.empty()) {
			throw std::runtime_error("Moving path file has no valid states: " + path);
		}
		return states;
	}

	void normalizeQuaternion(double q[4]) {
		const double n = std::sqrt(
			q[0] * q[0] + q[1] * q[1] + q[2] * q[2] + q[3] * q[3]);
		if (n <= kRangeEps) {
			throw std::runtime_error("Invalid zero-length moving target quaternion.");
		}
		for (int i = 0; i < 4; ++i) {
			q[i] /= n;
		}
	}

	void quaternionToMatrix(const double q[4], double R[3][3]) {
		const double w = q[0];
		const double x = q[1];
		const double y = q[2];
		const double z = q[3];

		R[0][0] = 1.0 - 2.0 * (y * y + z * z);
		R[0][1] = 2.0 * (x * y - z * w);
		R[0][2] = 2.0 * (x * z + y * w);

		R[1][0] = 2.0 * (x * y + z * w);
		R[1][1] = 1.0 - 2.0 * (x * x + z * z);
		R[1][2] = 2.0 * (y * z - x * w);

		R[2][0] = 2.0 * (x * z - y * w);
		R[2][1] = 2.0 * (y * z + x * w);
		R[2][2] = 1.0 - 2.0 * (x * x + y * y);
	}

	std::array<double, 3> normalizeVector(const std::array<double, 3>& v) {
		const double n = std::sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
		if (n <= kRangeEps) {
			throw std::runtime_error("Cannot normalize zero-length vector.");
		}
		return { v[0] / n, v[1] / n, v[2] / n };
	}

	std::array<double, 3> rotateGlobalToBody(
		const double R[3][3],
		const std::array<double, 3>& vGlobal)
	{
		return {
			R[0][0] * vGlobal[0] + R[1][0] * vGlobal[1] + R[2][0] * vGlobal[2],
			R[0][1] * vGlobal[0] + R[1][1] * vGlobal[1] + R[2][1] * vGlobal[2],
			R[0][2] * vGlobal[0] + R[1][2] * vGlobal[1] + R[2][2] * vGlobal[2]
		};
	}

	void vectorToThetaPhiDeg(
		const std::array<double, 3>& v,
		double& thetaDeg,
		double& phiDeg)
	{
		const std::array<double, 3> u = normalizeVector(v);
		const double z = std::max(-1.0, std::min(1.0, u[2]));
		thetaDeg = std::acos(z) * 180.0 / Pi;
		phiDeg = normalizePhiDegLocal(std::atan2(u[1], u[0]) * 180.0 / Pi);
		if (isThetaPoleLocal(thetaDeg)) {
			phiDeg = 0.0;
		}
	}

	MovingAngleSample computeMovingAngleSample(
		const MovingState& state,
		const std::array<double, 3>& radarPos)
	{
		const std::array<double, 3> rangeVec = {
			radarPos[0] - state.position[0],
			radarPos[1] - state.position[1],
			radarPos[2] - state.position[2]
		};
		const double range = std::sqrt(
			rangeVec[0] * rangeVec[0] +
			rangeVec[1] * rangeVec[1] +
			rangeVec[2] * rangeVec[2]);
		if (range <= kRangeEps) {
			throw std::runtime_error("Moving target range is too small.");
		}

		const std::array<double, 3> rGlobal = normalizeVector(rangeVec);
		double R[3][3];
		quaternionToMatrix(state.quat, R);
		const std::array<double, 3> rBody = rotateGlobalToBody(R, rGlobal);

		MovingAngleSample sample;
		sample.time_s = state.time_s;
		sample.pos[0] = state.position[0];
		sample.pos[1] = state.position[1];
		sample.pos[2] = state.position[2];
		sample.range = range;
		vectorToThetaPhiDeg(rGlobal, sample.thetaGlobalDeg, sample.phiGlobalDeg);
		vectorToThetaPhiDeg(rBody, sample.thetaBodyDeg, sample.phiBodyDeg);
		return sample;
	}

	std::vector<MovingAngleSample> computeMovingAngleSamples(
		const std::vector<MovingState>& states,
		const std::array<double, 3>& radarPos)
	{
		std::vector<MovingAngleSample> samples;
		samples.reserve(states.size());
		for (const auto& state : states) {
			samples.push_back(computeMovingAngleSample(state, radarPos));
		}
		return samples;
	}

	std::vector<MovingAngleSample> mergeMovingAngleSamples(
		const std::vector<MovingAngleSample>& samples,
		double mergeTolDeg)
	{
		if (mergeTolDeg <= 0.0) {
			throw std::runtime_error("Moving angle merge tolerance must be positive.");
		}

		std::vector<MovingAngleSample> uniqueSamples;
		for (auto sample : samples) {
			if (isThetaPoleLocal(sample.thetaBodyDeg)) {
				sample.phiBodyDeg = 0.0;
			}

			bool exists = false;
			for (const auto& kept : uniqueSamples) {
				if (std::abs(sample.thetaBodyDeg - kept.thetaBodyDeg) <= mergeTolDeg &&
					deltaPhiDeg(sample.phiBodyDeg, kept.phiBodyDeg) <= mergeTolDeg) {
					exists = true;
					break;
				}
			}
			if (!exists) {
				uniqueSamples.push_back(sample);
			}
		}
		return uniqueSamples;
	}
}
