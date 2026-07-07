// OCTree.cpp
#include "OCTree.h"
#include "PathUtils.h"
#include "OverloadAlgo.h"

#include <algorithm>
#include <functional>
#include <fstream>
#include <iostream>
#include <limits>
#include <unordered_set>

OCTree::OCTree(const std::vector<RWGBase>& rwgs, double wavelength, const std::string& nasPath)
	: nasPath(nasPath)
{
	computeBoundingBox(rwgs, wavelength);

	buildLeafNodes(rwgs);

	buildUpperLevels();

	assignNodeNumbers();

	buildNearNeighbors();

	buildDistantRelatives();

	verifyNeighborsAndDR();

	/*exportOctreeToFile();
	exportNodeToFile();*/
}

// 导出远亲近邻信息
void OCTree::exportOctreeToFile() {
	std::string parentPath = getParentPath(nasPath);
	std::string baseName = getBaseName(nasPath);
	std::string currentTime = getCurrentTimeString();
	std::ofstream infoFile(parentPath + baseName + "_OCTree_near_distant_" + currentTime + "_.txt");
	if (!infoFile.is_open()) {
		std::cerr << "Error: Cannot open file: " << nasPath << " (OCTree.cpp)" << std::endl;
		return;
	}

	infoFile << std::fixed << std::setprecision(12);
	infoFile << "节点信息\n";
	for (int i = 1; i <= maxLevel_; ++i) {
		if (i == 1) {
			for (int j = 0; j < octreeNodes_[i].size(); ++j) {
				Node* node = octreeNodes_[i][j];
				infoFile << "第 " << i << " 层: " << octreeNodes_[i].size() << "\n";
				infoFile << i << "\t" << node->nodeNumber << "\t" << " 中心点: ( "
					<< node->center.x << " , " << node->center.y << " , " << node->center.z
					<< " )， RWG数量为: " << node->rwgIndices.size() << "\n";
				for (int k = 0; k < node->nearNeighbors.size(); ++k) {
					infoFile << "第 " << k + 1 << " 个，近邻序号：" << node->nearNeighbors[k]->nodeNumber << "\n";
				}
				infoFile << "\n";
			}
			infoFile << "\n\n";
		}
		else {
			for (int j = 0; j < octreeNodes_[i].size(); ++j) {
				Node* node = octreeNodes_[i][j];
				infoFile << "第 " << i << " 层: " << octreeNodes_[i].size() << "\n";
				infoFile << i << "\t" << node->nodeNumber << "\t" << " 中心点: ( "
					<< node->center.x << " , " << node->center.y << " , " << node->center.z
					<< " )， RWG数量为: " << node->rwgIndices.size() << "\n";
				for (int k = 0; k < node->nearNeighbors.size(); ++k) {
					infoFile << "第 " << k + 1 << " 个，近邻序号：" << node->nearNeighbors[k]->nodeNumber << "\n";
				}
				infoFile << "\n";
				for (int l = 0; l < node->distantRelatives.size(); ++l) {
					infoFile << "第 " << l + 1 << " 个，远亲序号：" << node->distantRelatives[l]->nodeNumber << " 向量为: ( "
						<< node->center.x - node->distantRelatives[l]->center.x << " , "
						<< node->center.y - node->distantRelatives[l]->center.y << " , "
						<< node->center.z - node->distantRelatives[l]->center.z << " )" << "\n";
					infoFile << "drid: " << node->dRVecId[l] << "\n";
					infoFile << "drVec: ( " << octreeNodesDRvec_[i - 2][node->dRVecId[l]].x << " , "
						<< octreeNodesDRvec_[i - 2][node->dRVecId[l]].y << " , "
						<< octreeNodesDRvec_[i - 2][node->dRVecId[l]].z << " )\n";
					infoFile << "\n";
				}
				infoFile << "\n\n";
			}
		}
		infoFile << "\n\n";
	}
	infoFile.close();
	std::cout << "Info: Exported to ==> " << parentPath + baseName + "_OCTree_near_distant_" + currentTime + "_.txt (OCTree.cpp)" << std::endl;
}

void OCTree::exportNodeToFile() {
	std::string parentPath = getParentPath(nasPath);
	std::string baseName = getBaseName(nasPath);
	std::string currentTime = getCurrentTimeString();
	std::ofstream infoFile(parentPath + baseName + "_OCTree_nodes_" + currentTime + "_.txt");
	if (!infoFile.is_open()) {
		std::cerr << "Error: Cannot open file: " << nasPath << " (OCTree.cpp)" << std::endl;
		return;
	}

	// ===== 2. 输出 =====
	for (int lvl = 0; lvl < octreeNodes_.size(); ++lvl) {
		const auto& level = octreeNodes_[lvl];

		infoFile << "========== Level " << lvl
			<< " | Node count = " << level.size()
			<< " ==========\n";

		for (auto node : level) {
			infoFile << "NodeID: " << node->nodeNumber << "\n";

			infoFile << "Center: ( "
				<< node->center.x << ", "
				<< node->center.y << ", "
				<< node->center.z << " )\n";

			infoFile << "Morton: " << node->mortonCode << "\n";
			infoFile << "Level: " << node->level << "\n";
			infoFile << "Cube Length: " << node->cubeLength << "\n";
			infoFile << "Parent Morton: " << node->parentCode << "\n";

			infoFile << "Father ID: "
				<< (node->realFather ? node->realFather->nodeNumber : -1)
				<< "\n";

			infoFile << "Children count: " << node->realSons.size() << "\n";

			infoFile << "RWG count: " << node->rwgIndices.size() << "\n";

			infoFile << "Near count: " << node->nearNeighbors.size() << "\n";
			infoFile << "Distant count: " << node->distantRelatives.size() << "\n";

			infoFile << "vecToFather: ( "
				<< node->vecToFather.x << ", "
				<< node->vecToFather.y << ", "
				<< node->vecToFather.z << " )\n";

			// ===== 邻居ID =====
			infoFile << "Near IDs: ";
			for (auto n : node->nearNeighbors)
				infoFile << n->nodeNumber << " ";
			infoFile << "\n";

			infoFile << "Distant IDs: ";
			for (auto n : node->distantRelatives)
				infoFile << n->nodeNumber << " ";
			infoFile << "\n";

			infoFile << "----------------------------------\n";
		}

		infoFile << "\n\n";
	}

	infoFile.close();
	std::cout << "Info: Exported to ==> " << parentPath + baseName + "_OCTree_nodes_" + currentTime + "_.txt (OCTree.cpp)" << std::endl;
}

void OCTree::computeBoundingBox(const std::vector<RWGBase>& rwgs_, const double wavelength) {
	// 01. 初始化
	globalBBoxMin_ = {
		std::numeric_limits<double>::max(),
		std::numeric_limits<double>::max(),
		std::numeric_limits<double>::max()
	};

	globalBBoxMax_ = {
		std::numeric_limits<double>::lowest(),
		std::numeric_limits<double>::lowest(),
		std::numeric_limits<double>::lowest()
	};

	// 02. 计算包围盒
	for (const auto& rwg : rwgs_) {
		const Point& p = rwg.edgeCenterCenter;

		globalBBoxMin_.x = std::min(globalBBoxMin_.x, p.x);
		globalBBoxMin_.y = std::min(globalBBoxMin_.y, p.y);
		globalBBoxMin_.z = std::min(globalBBoxMin_.z, p.z);

		globalBBoxMax_.x = std::max(globalBBoxMax_.x, p.x);
		globalBBoxMax_.y = std::max(globalBBoxMax_.y, p.y);
		globalBBoxMax_.z = std::max(globalBBoxMax_.z, p.z);
	}

	double dx = globalBBoxMax_.x - globalBBoxMin_.x;
	double dy = globalBBoxMax_.y - globalBBoxMin_.y;
	double dz = globalBBoxMax_.z - globalBBoxMin_.z;

	double rootLen = std::max(dx, std::max(dy, dz)) * 1.001;

	// 03. 计算leaf尺寸
	double leafTarget = leafCoeff_ * wavelength;

	double ratio = rootLen / leafTarget;
	if (ratio < 1.0) ratio = 1.0;

	maxLevel_ = static_cast<int>(std::ceil(std::log2(ratio)));
	maxLevel_ = std::max(2, std::min(maxLevel_, 15));

	rootLen = leafTarget * (1 << maxLevel_);

	// 04. 计算中心
	double cx = 0.5 * (globalBBoxMin_.x + globalBBoxMax_.x);
	double cy = 0.5 * (globalBBoxMin_.y + globalBBoxMax_.y);
	double cz = 0.5 * (globalBBoxMin_.z + globalBBoxMax_.z);

	// 05. 重新对齐bbox(变成立方体)
	globalBBoxMin_ = {
		cx - 0.5 * rootLen,
		cy - 0.5 * rootLen,
		cz - 0.5 * rootLen
	};

	globalBBoxMax_ = {
		cx + 0.5 * rootLen,
		cy + 0.5 * rootLen,
		cz + 0.5 * rootLen
	};

	leafLen_ = rootLen / (1 << maxLevel_);

	// 06. 构造root节点
	rootNode_ = Node(
		0,                      // mortonCode
		0,                      // level
		{ cx, cy, cz },         // center
		rootLen,                // cubeLength
		0,                      // parentCode
		{ 0,0,0 }               // vecToFather
	);

	// 09. 输出信息
	std::cout << "==================================" << std::endl;
	std::cout << "BBox size: " << rootLen << " m" << std::endl;
	std::cout << "lambda: " << wavelength << " m" << std::endl;
	std::cout << "leafCoeff: " << leafCoeff_ << std::endl;
	std::cout << "maxLevel: " << maxLevel_ << std::endl;
	std::cout << "leaf size (λ): " << leafLen_ / wavelength << std::endl;
}

void OCTree::buildLeafNodes(const std::vector<RWGBase>& rwg_) {
	std::unordered_map<uint64_t, std::unique_ptr<Node>> leafMap;

	for (const auto& rwg : rwg_) {
		const auto& p = rwg.edgeCenterCenter;

		int ix = static_cast<int>((p.x - globalBBoxMin_.x) / leafLen_);
		int iy = static_cast<int>((p.y - globalBBoxMin_.y) / leafLen_);
		int iz = static_cast<int>((p.z - globalBBoxMin_.z) / leafLen_);

		// 防止越界
		int maxIdx = (1 << maxLevel_) - 1;
		ix = std::max(0, std::min(ix, maxIdx));
		iy = std::max(0, std::min(iy, maxIdx));
		iz = std::max(0, std::min(iz, maxIdx));

		// ===== 2. Morton编码 =====
		uint64_t morton = morton3D(ix, iy, iz);

		// ===== 3. 插入或查找节点 =====
		auto& nodePtr = leafMap[morton];

		if (!nodePtr) {
			double cx = globalBBoxMin_.x + (ix + 0.5) * leafLen_;
			double cy = globalBBoxMin_.y + (iy + 0.5) * leafLen_;
			double cz = globalBBoxMin_.z + (iz + 0.5) * leafLen_;

			nodePoint fatherCenter{
				globalBBoxMin_.x + ((ix >> 1) + 0.5) * (2 * leafLen_),
				globalBBoxMin_.y + ((iy >> 1) + 0.5) * (2 * leafLen_),
				globalBBoxMin_.z + ((iz >> 1) + 0.5) * (2 * leafLen_)
			};

			nodePoint vecToFather{
				fatherCenter.x - cx,
				fatherCenter.y - cy,
				fatherCenter.z - cz
			};

			nodePtr = std::make_unique<Node>(
				morton,
				maxLevel_,
				nodePoint{ cx, cy, cz },
				leafLen_,
				morton >> 3,
				vecToFather
			);
		}

		nodePtr->rwgIndices.push_back(const_cast<RWGBase*>(&rwg));
	}
	// ===== 5. 存入 octreeNodes_ =====
	octreeNodes_.resize(maxLevel_ + 1);

	for (auto& [code, node] : leafMap) {
		Node* rawPtr = node.get();
		octreeNodes_[maxLevel_].push_back(node.get());
		nodePool_.push_back(std::move(node));
	}

	std::cout << "Leaf nodes count = "
		<< octreeNodes_[maxLevel_].size()
		<< std::endl;
}

void OCTree::buildUpperLevels() {
	// 确保每一层vector存在
	octreeNodes_.resize(maxLevel_ + 1);

	for (int level = maxLevel_; level > 0; --level) {
		std::unordered_map<uint64_t, Node*> parentMap;

		// ===== 1. 遍历当前层（child层）=====
		for (Node* child : octreeNodes_[level]) {
			uint64_t parentCode = child->mortonCode >> 3;
			int childIdx = child->mortonCode & 0b111;

			Node*& parent = parentMap[parentCode];

			// ===== 2. 如果父节点不存在 → 创建 =====
			if (!parent) {
				double parentLen = child->cubeLength * 2.0;
				double offset = child->cubeLength * 0.5;

				nodePoint parentCenter = child->center;

				// 关键：反推父节点中心
				parentCenter.x += (childIdx & 1) ? -offset : offset;
				parentCenter.y += (childIdx & 2) ? -offset : offset;
				parentCenter.z += (childIdx & 4) ? -offset : offset;

				auto newNode = std::make_unique<Node>(
					parentCode,
					level - 1,
					parentCenter,
					parentLen,
					parentCode >> 3,
					nodePoint{ 0,0,0 }
				);

				parent = newNode.get();
				nodePool_.push_back(std::move(newNode));
			}

			// ===== 3. 建立父子关系 =====
			child->realFather = parent;
			child->vecToFather = {
				parent->center.x - child->center.x,
				parent->center.y - child->center.y,
				parent->center.z - child->center.z
			};

			parent->realSons.push_back(child);

			// ===== 4. 聚合RWG =====
			parent->rwgIndices.insert(
				parent->rwgIndices.end(),
				child->rwgIndices.begin(),
				child->rwgIndices.end()
			);
		}

		// ===== 5. 存入上一层 =====
		octreeNodes_[level - 1].reserve(parentMap.size());
		for (auto& [code, nodePtr] : parentMap) {
			octreeNodes_[level - 1].push_back(nodePtr);
		}

		std::cout << "Level " << level - 1 << " node count = "
			<< octreeNodes_[level - 1].size() << std::endl;
	}
}

void OCTree::assignNodeNumbers() {
	for (int level = 0; level <= maxLevel_; ++level) {
		for (int num = 0; num < octreeNodes_[level].size(); ++num) {
			octreeNodes_[level][num]->nodeNumber = num;
		}
	}
}

void OCTree::buildNearNeighbors() {
	for (int level = 0; level <= maxLevel_; ++level) {
		// 建一个 Morton → Node* 映射（加速查找）
		std::unordered_map<uint64_t, Node*> nodeMap;

		for (Node* node : octreeNodes_[level]) {
			nodeMap[node->mortonCode] = node;
			node->nearNeighbors.clear();
		}

		int gridSize = 1 << level;

		for (Node* node : octreeNodes_[level]) {
			int ix, iy, iz;
			decodeMorton3D(node->mortonCode, ix, iy, iz);

			for (int dx = -1; dx <= 1; ++dx)
				for (int dy = -1; dy <= 1; ++dy)
					for (int dz = -1; dz <= 1; ++dz) {
						if (dx == 0 && dy == 0 && dz == 0)continue;

						int nx = ix + dx;
						int ny = iy + dy;
						int nz = iz + dz;

						if (nx < 0 || ny < 0 || nz < 0 ||
							nx >= gridSize || ny >= gridSize || nz >= gridSize)
							continue;

						uint64_t code = morton3D(nx, ny, nz);

						auto it = nodeMap.find(code);
						if (it != nodeMap.end()) {
							node->nearNeighbors.push_back(it->second);
						}
					}
		}

		std::cout << "Level " << level
			<< " near neighbors built." << std::endl;
	}
}

void OCTree::buildDistantRelatives() {
	for (int level = 2; level <= maxLevel_; ++level) {
		std::unordered_map<uint64_t, Node*> nodeMap;

		for (Node* node : octreeNodes_[level]) {
			nodeMap[node->mortonCode] = node;
		}

		for (Node* node : octreeNodes_[level]) {
			Node* parent = node->realFather;
			if (!parent)continue;

			std::unordered_set<Node*> drSet;

			int px, py, pz;
			decodeMorton3D(parent->mortonCode, px, py, pz);

			int gridSize = 1 << (level - 1);

			for (int dx = -1; dx <= 1; ++dx)
				for (int dy = -1; dy <= 1; ++dy)
					for (int dz = -1; dz <= 1; ++dz) {
						int nx = px + dx;
						int ny = py + dy;
						int nz = pz + dz;

						if (nx < 0 || ny < 0 || nz < 0 ||
							nx >= gridSize || ny >= gridSize || nz >= gridSize)
							continue;

						// 👉 遍历这个 parent cell 的 8 个 child
						for (int cx = 0; cx < 2; ++cx)
							for (int cy = 0; cy < 2; ++cy)
								for (int cz = 0; cz < 2; ++cz)
								{
									int child_x = nx * 2 + cx;
									int child_y = ny * 2 + cy;
									int child_z = nz * 2 + cz;

									uint64_t childCode = morton3D(child_x, child_y, child_z);

									auto it = nodeMap.find(childCode);  // 当前层查找
									if (it == nodeMap.end()) continue;

									Node* child = it->second;

									if (child == node) continue;

									if (!isNear(node, child)) {
										drSet.insert(child);
									}
								}
					}
			node->distantRelatives.assign(drSet.begin(), drSet.end());
		}
		std::cout << "Level " << level
			<< " Distant Relatives built." << std::endl;
	}

	getDRVec();
}

bool OCTree::isNear(Node* a, Node* b) {
	double dist = std::max({
		std::fabs(a->center.x - b->center.x),
		std::fabs(a->center.y - b->center.y),
		std::fabs(a->center.z - b->center.z)
		});

	return dist <= (1.74 * a->cubeLength);
}

void OCTree::getDRVec() {
	// 提前校验，避免越界
	if (maxLevel_ < 2) {
		octreeNodesDRvec_.clear();
		return;
	}
	octreeNodesDRvec_.resize(maxLevel_ - 1);

	for (int i = 2; i <= maxLevel_; ++i) {
		auto& drVecLevel = octreeNodesDRvec_[i - 2];
		const auto& levelNodes = octreeNodes_[i];
		// 哈希表
		std::unordered_map<VecKey, int, VecKeyHash> vecMap;

		for (Node* node : levelNodes) {
			if (node->distantRelatives.empty())continue;
			node->dRVecId.clear();

			for (Node* dr : node->distantRelatives) {
				// ===== 1. 转整数位移 =====
				int dx = int(std::round((node->center.x - dr->center.x) / node->cubeLength));
				int dy = int(std::round((node->center.y - dr->center.y) / node->cubeLength));
				int dz = int(std::round((node->center.z - dr->center.z) / node->cubeLength));

				VecKey key{ dx, dy, dz };

				// ===== 2. 查找或插入 =====
				auto it = vecMap.find(key);

				int vecIdx;

				if (it == vecMap.end()) {
					vecIdx = drVecLevel.size();

					vecMap[key] = vecIdx;

					// 同时存真实向量（用于后续计算）
					drVecLevel.push_back({
						dx * node->cubeLength,
						dy * node->cubeLength,
						dz * node->cubeLength
						});
				}
				else {
					vecIdx = it->second;
				}

				node->dRVecId.push_back(vecIdx);
			}
		}
	}
}

void OCTree::verifyNeighborsAndDR() {
	std::cout << "\n========== Verifying Near & DR ==========\n";

	for (int level = 2; level <= maxLevel_; ++level) {
		const auto& nodes = octreeNodes_[level];

		int nearError = 0;
		int drError = 0;

		for (Node* a : nodes) {
			std::unordered_set<Node*> nearTruth;
			std::unordered_set<Node*> drTruth;

			// ===== 1. 几何 brute force =====
			for (Node* b : nodes) {
				if (a == b) continue;

				bool nearGeom = cubeNear(a, b);

				if (nearGeom) {
					nearTruth.insert(b);
				}
				else {
					// 判断是否为远亲
					Node* pa = a->realFather;
					Node* pb = b->realFather;

					if (pa && pb && cubeNear(pa, pb)) {
						drTruth.insert(b);
					}
				}
			}

			// ===== 2. 验证 near =====
			std::unordered_set<Node*> nearBuilt(
				a->nearNeighbors.begin(),
				a->nearNeighbors.end()
			);

			for (auto* n : nearTruth) {
				if (!nearBuilt.count(n)) {
					nearError++;
				}
			}

			for (auto* n : nearBuilt) {
				if (!nearTruth.count(n)) {
					nearError++;
				}
			}

			// ===== 3. 验证 DR =====
			std::unordered_set<Node*> drBuilt(
				a->distantRelatives.begin(),
				a->distantRelatives.end()
			);

			for (auto* d : drTruth) {
				if (!drBuilt.count(d)) {
					drError++;
				}
			}

			for (auto* d : drBuilt) {
				if (!drTruth.count(d)) {
					drError++;
				}
			}
		}

		std::cout << "Level " << level
			<< "  Near Error = " << nearError
			<< "  DR Error = " << drError
			<< std::endl;
	}

	std::cout << "=========================================\n";
}

bool OCTree::cubeNear(Node* a, Node* b) {
	double dx = std::fabs(a->center.x - b->center.x);
	double dy = std::fabs(a->center.y - b->center.y);
	double dz = std::fabs(a->center.z - b->center.z);

	double len = 1.01 * a->cubeLength;

	return (dx <= len && dy <= len && dz <= len);
}