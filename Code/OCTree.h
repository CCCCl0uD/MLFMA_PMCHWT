// OCTree.h
#pragma once
#ifndef OCTREE_H
#define OCTREE_H

#include"MyStruct.h"

#include <vector>
#include <cstdint>
#include <memory>
#include <string>

class OCTree {
public:
	using nodePoint = ::nodePoint;

	static inline double snapZero(double x) {
		return(std::abs(x) < 1e-6) ? 0.0 : x;
	}

	static inline bool isVecEqual(const OCTree::nodePoint& a, const OCTree::nodePoint& b) {
		return std::abs(a.x - b.x) < 1e-6 &&
			std::abs(a.y - b.y) < 1e-6 &&
			std::abs(a.z - b.z) < 1e-6;
	}

	struct VecKey {
		int dx, dy, dz;
		bool operator==(const VecKey& other) const {
			return dx == other.dx &&
				dy == other.dy &&
				dz == other.dz;
		}
	};

	struct VecKeyHash {
		std::size_t operator()(const VecKey& k) const {
			return((std::hash<int>()(k.dx) ^ (std::hash<int>()(k.dy) << 1)) >> 1)
				^ (std::hash<int>()(k.dz) << 1);
		}
	};

	// 节点结构体
	struct Node {
		int nodeNumber = -1;                // 节点序号
		int level = 0;                      // 节点层级
		double cubeLength = 0.0;            // 该盒子边长
		nodePoint center;                   // 节点中心坐标

		uint64_t mortonCode = 0;            // 该盒子的Morton码
		uint64_t parentCode = 0;            // 直系父亲的Morton码

		std::vector<RWGBase*> rwgIndices;   // 存储的RWG索引

		Node* realFather = nullptr;         // 直系父亲的指针
		nodePoint vecToFather;              // 该盒子到父亲的向量
		std::vector<Node*> realSons;        // 直系子族的指针

		std::vector<Node*> nearNeighbors;   // 近邻节点的指针
		std::vector<Node*> distantRelatives;// 远亲节点的指针

		std::vector<int> dRVecId;           // 远亲节点的向量序号

		Node() = default;

		Node(uint64_t mortonCode,
			int level,
			nodePoint center,
			double cubeLength,
			uint64_t parentCode,
			nodePoint vecToFather)
			: mortonCode(mortonCode),
			level(level),
			center(center),
			cubeLength(cubeLength),
			parentCode(parentCode),
			vecToFather(vecToFather) {
		}
	};

	OCTree(const std::vector<RWGBase>& rwgs, double wavelength, const std::string& nasPath);
	const std::string nasPath;
	int maxLevel_;
	std::vector<std::vector<Node*>> octreeNodes_;
	std::vector<std::vector<nodePoint>> octreeNodesDRvec_; // 远亲节点的向量 (2, 3, ...., maxLevel_ -1, maxLevel_)

private:

	void exportOctreeToFile();
	void exportNodeToFile();
	void computeBoundingBox(const std::vector<RWGBase>& rwgs_, const double wavelength);
	void buildLeafNodes(const std::vector<RWGBase>& rwg_);
	void buildUpperLevels();
	void assignNodeNumbers();
	void buildNearNeighbors();
	void buildDistantRelatives();
	bool isNear(Node* a, Node* b);
	void getDRVec();
	void verifyNeighborsAndDR();
	bool cubeNear(Node* a, Node* b);

	nodePoint globalBBoxMin_;
	nodePoint globalBBoxMax_;
	double leafLen_;

	Node rootNode_;

	std::vector<std::unique_ptr<Node>> nodePool_;

	inline uint64_t expandBits(uint32_t x) {
		uint64_t v = x;
		v = (v | (v << 32)) & 0x1f00000000ffff;
		v = (v | (v << 16)) & 0x1f0000ff0000ff;
		v = (v | (v << 8)) & 0x100f00f00f00f00f;
		v = (v | (v << 4)) & 0x10c30c30c30c30c3;
		v = (v | (v << 2)) & 0x1249249249249249;
		return v;
	}

	inline uint64_t morton3D(uint32_t x, uint32_t y, uint32_t z) {
		return (expandBits(x)) |
			(expandBits(y) << 1) |
			(expandBits(z) << 2);
	}

	inline uint32_t compactBits(uint64_t x) {
		x &= 0x1249249249249249ULL;              // 只保留 x 位
		x = (x ^ (x >> 2)) & 0x10c30c30c30c30c3ULL;
		x = (x ^ (x >> 4)) & 0x100f00f00f00f00fULL;
		x = (x ^ (x >> 8)) & 0x1f0000ff0000ffULL;
		x = (x ^ (x >> 16)) & 0x1f00000000ffffULL;
		x = (x ^ (x >> 32)) & 0x00000000001fffffULL;
		return static_cast<uint32_t>(x);
	}

	// Morton解码
	inline void decodeMorton3D(uint64_t code, int& ix, int& iy, int& iz) {
		ix = compactBits(code);
		iy = compactBits(code >> 1);
		iz = compactBits(code >> 2);
	}
};

#endif // OCTREE_H