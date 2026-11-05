// main.cpp
#include "BatchRunner.h"
#include "MoM.h"
#include "FMM.h"
#include "MLFMM.h"
#include "RWGGenerator.h"
#include "RunLogger.h"

//\f:Times New Roman(f = 0.8GHz)
//\f:Times New Roman(HH polarizations)
//\g(q)\ - (\f:Times New Roman(i)) = \g(q)\ - (\f:Times New Roman(s)) = 90\ + (�?, \g(f)\ - (\f:Times New Roman(i)) = 0\ + (�?
//\g(e)\ - (\f:Times New Roman(r)) = (4, -0.001)

int main(int argc, char* argv[]) {
	// ------- Batch -------
	if (argc >= 2) {
		std::string arg1 = argv[1];
		if (arg1 == "batch") {
			if (argc < 3) {
				std::cerr << "Usage: MLFMA.exe batch <config_file>\n";
				return 1;
			}
			try {
				BatchRunner::run(argv[2]);
			}
			catch (const std::exception& e) {
				std::cerr << "Batch failed: " << e.what() << "\n";
				return 1;
			}
			return 0;
		}
		else if (arg1 == "--help" || arg1 == "-h") {
			std::cout << "Usage:\n"
				<< "  MLFMA.exe              (single mode)\n"
				<< "  MLFMA.exe batch <cfg>  (batch mode)\n";
			return 0;
		}
		else {
			std::cerr << "Unknown argument: " << arg1 << "\n";
			std::cerr << "Usage: MLFMA.exe [batch <config_file>]\n";
			return 1;
		}
	}

	// ------- Single -------
	std::vector<Point> points;
	std::vector<FaceElement> faces;
	const std::string inputFile = "D:\\MyCode\\PMCHWT_MLFMA\\DATA\\Muti\\4_Sphere\\4_Sphere_1e9.nas";
	RCSExportConfig cfg(inputFile);
	cfg.movingPathFile = "D:\\MyCode\\PMCHWT_MLFMA\\DATA\\Cube_0d5_0f20\\moving_path_Cube_1e9_z.txt";
	const double freq = 1.0e9, E0 = 1.0;
	cfg.incTheta = 90.0;
	cfg.incPhi = 0.0;
	cfg.scaThetaStart = 90.0;
	cfg.scaThetaEnd = 90.0;
	cfg.scaPhiStart = 0.0;
	cfg.scaPhiEnd = 360.0;
	cfg.scaStep = 1.0;
	cfg.nPoints = 7; // 1, 3, 4, 6, 7, 12, 13, 16

	cfg.selectAlgo = selectAlgorithm::MLFMA;
	cfg.integralEqu = IntegralEquation::CFIE;
	cfg.monoDual = MonoDual::Dual;
	cfg.polWave = Polarization::V;
	const int selectMatrixSolver = 0; // 0 = GMRES, 1 = CGS

	cfg.useFSBMono = false;
	cfg.fsbRho = 1.2;
	cfg.exportFSBComplexField = false;
	cfg.exportFSBSamples = false;
	cfg.radarX = 0.5;
	cfg.radarY = 0.5;
	cfg.radarZ = 0.7071067812;
	cfg.useFSBMLFMAMoving = false;
	cfg.movingAngleMarginDeg = 2.0;
	cfg.movingAngleMergeTolDeg = 0.25;
	cfg.exportMovingAngles = false;
	cfg.exportMovingComplexField = false;
	cfg.exportMovingSamples = false;
	cfg.useAccurateMLFMAMoving = false;
	cfg.useACARHSMoving = false;
	cfg.acaRHSTol = 1.0e-3;
	cfg.acaRHSMaxRank = 0;
	cfg.acaRHSPivotEps = 1.0e-14;
	cfg.acaRHSCheckCount = 10;
	cfg.exportACARHSCheck = false;
	cfg.exportACARHSFactors = false;

	if (cfg.integralEqu == IntegralEquation::PMCHWT) {
		cfg.material = materialType::Die;
		cfg.epsilonR = { 1.0, -0.01 };
		cfg.muR = { 1.0, 0.0 };
	}

	omp_set_dynamic(0);
	omp_set_num_threads(12);
	/************************************************************************/
	std::ios::sync_with_stdio(false);
	std::cin.tie(nullptr);
	/************************************************************************/
	int debugInfoToFile = 0;// 0=>off�?=>on
	/************************************************************************/
	if (selectMatrixSolver != 0 && selectMatrixSolver != 1) {
		std::cerr << "Error: Invalid selectMatrixSolver value: "
			<< selectMatrixSolver
			<< " (0 = GMRES, 1 = CGS)" << std::endl;
		return 3;
	}

	cfg.validate();
	EMSource wave(freq, cfg.polWave == Polarization::V ? 0.0 : 90.0);
	if (cfg.material == materialType::Die) {
		wave.initDielectric(cfg.epsilonR, cfg.muR);
	}
	ScopedRunLogger runLogger(cfg.buildFilePath("_RUN_LOG", "txt"));
	if (runLogger.enabled()) {
		std::cout << "Run log file: " << runLogger.path() << std::endl;
	}

	gaussPoints gausspoint(cfg.nPoints);

	MeshReadReport meshReport;
	if (readNasData(inputFile, points, faces, meshReport) != 0) {
		std::cerr << "Error: Failed to read Nastran file (main.cpp)" << std::endl;
		return 1;
	}

	std::vector<RWGBase> rwgBases;
	RWGOptions rwgOptions;
	rwgOptions.requireClosedSurface = cfg.integralEqu != IntegralEquation::EFIE;
	RWGGenerationReport rwgReport;
	if (!RWGGenerator::GenerateRWG(faces, points, rwgBases, rwgOptions, rwgReport)) {
		std::cerr << "Error: RWG basis generation failed (main.cpp)" << std::endl;
		return 2;
	}

	if (debugInfoToFile == 1) {
		if (!writeOutputFiles(inputFile, faces, points)) {
			std::cerr << "Error: Cannot write face element data to file (main.cpp)" << std::endl;
			return 1;
		}

		if (!RWGGenerator::WriteRWGToFile(rwgBases, inputFile)) {
			std::cerr << "Error: Cannot write RWG basis data to file (main.cpp)" << std::endl;
			return 2;
		}
	}

	if (cfg.selectAlgo == selectAlgorithm::MoM) {
		auto startMoM = std::chrono::high_resolution_clock::now();
		MoM mom(cfg, selectMatrixSolver, rwgBases, gausspoint, E0, wave);
		auto endMoM = std::chrono::high_resolution_clock::now();
		auto durationMoM = std::chrono::duration_cast<std::chrono::microseconds>(endMoM - startMoM);
		std::cout << "<<<< Info: MoM done, elapsed: " << durationMoM.count() / 1000.0 << " ms >>>> (main.cpp)\n" << std::endl;
		return 0;
	}

	auto startOCTree = std::chrono::high_resolution_clock::now();
	OCTree octree(rwgBases, wave.wavelength(), inputFile);
	auto endOCTree = std::chrono::high_resolution_clock::now();
	auto durationOCTree = std::chrono::duration_cast<std::chrono::microseconds>(endOCTree - startOCTree);
	std::cout << "<<<< Info: OCTree built, elapsed: " << durationOCTree.count() / 1000.0 << " ms >>>> (main.cpp)\n" << std::endl;

	if (cfg.selectAlgo == selectAlgorithm::FMM) {
		auto startFMM = std::chrono::high_resolution_clock::now();
		FMM fmm(cfg, selectMatrixSolver,
			octree.octreeNodes_, rwgBases, octree.maxLevel_, gausspoint,
			E0, wave);
		auto endFMM = std::chrono::high_resolution_clock::now();
		auto durationFMM = std::chrono::duration_cast<std::chrono::microseconds>(endFMM - startFMM);
		std::cout << "<<<< Info: FMM done, elapsed: " << durationFMM.count() / 1000.0 << " ms >>>> (main.cpp)\n" << std::endl;
	}
	else if (cfg.selectAlgo == selectAlgorithm::MLFMA) {
		auto startMLFMA = std::chrono::high_resolution_clock::now();
		MLFMM mlfmm(cfg, selectMatrixSolver,
			octree.octreeNodes_, octree.octreeNodesDRvec_, rwgBases, octree.maxLevel_, gausspoint,
			E0, wave);
		auto endMLFMA = std::chrono::high_resolution_clock::now();
		auto durationMLFMA = std::chrono::duration_cast<std::chrono::microseconds>(endMLFMA - startMLFMA);
		std::cout << "<<<< Info: MLFMA done, elapsed: " << durationMLFMA.count() / 1000.0 << " ms >>>> (main.cpp)\n" << std::endl;
	}

	return 0;
}