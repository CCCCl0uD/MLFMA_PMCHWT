// RWGGenerator.h
#pragma once
#ifndef RWGGENERATOR_H
#define RWGGENERATOR_H
#include "NASReader.h"
#include <stdexcept>
#include <string>
#include <vector>

class RWGGenerator {
public:
	// 生成RWG基函数数据
	static bool GenerateRWG(const std::vector<FaceElement>& faces,
		const std::vector<Point>& points,
		std::vector<RWGBase>& result);

	// 将RWG基函数数据写入文件
	static bool WriteRWGToFile(const std::vector<RWGBase>& bases,
		const std::string& nasPath);

private:
	// 边描述结构体
	struct Edge {
		int v[2];  // 排序后的顶点ID

		Edge(int a, int b) {
			v[0] = (a < b) ? a : b;
			v[1] = (a < b) ? b : a;
		}

		bool operator<(const Edge& other) const {
			return (v[0] < other.v[0]) || (v[0] == other.v[0] && v[1] < other.v[1]);
		}
	};

	// 带面元索引的边信息
	struct EdgeInfo {
		Edge edge;
		size_t faceIndex;

		EdgeInfo(int a, int b, size_t idx)
			: edge(a, b), faceIndex(idx) {
		}
	};
};

// 辅助函数声明
namespace PathUtil {
	std::string GetParentPath(const std::string& fullPath);
	std::string GetBaseName(const std::string& fullPath);
}

#endif //RWGGENERATOR_H