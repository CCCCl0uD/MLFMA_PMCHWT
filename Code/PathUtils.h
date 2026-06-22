#pragma once
#pragma once
#ifndef PATHUTILS_H
#define PATHUTILS_H

#include <chrono>
#include <ctime>
#include <iomanip>
#include <sstream>
#include <string>

inline std::string getCurrentTimeString() {
	auto now = std::chrono::system_clock::now();
	std::time_t now_time = std::chrono::system_clock::to_time_t(now);
	std::tm local_tm;
	localtime_s(&local_tm, &now_time);
	std::ostringstream oss;
	oss << std::put_time(&local_tm, "%Y%m%d_%H%M%S");
	return oss.str();
}

inline std::string getParentPath(const std::string& filePath) {
	size_t pos = filePath.find_last_of("\\/");
	if (pos != std::string::npos) {
		return filePath.substr(0, pos);
	}
	return ".";
}

inline std::string getBaseName(const std::string& filePath) {
	std::string name = filePath;
	size_t pos1 = name.find_last_of("\\/");
	if (pos1 != std::string::npos) {
		name = name.substr(pos1 + 1);
	}
	size_t pos2 = name.find_last_of('.');
	if (pos2 != std::string::npos) {
		name = name.substr(0, pos2);
	}
	return name;
}

#endif // PATHUTILS_H