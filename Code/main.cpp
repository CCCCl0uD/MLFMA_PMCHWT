// main.cpp
#include "BatchRunner.h"
#include "MoM.h"
#include "FMM.h"
#include "MLFMM.h"
#include "RWGGenerator.h"
#include "RunLogger.h"

RCSExportConfig createRCSExportConfig(const std::string& inputFile,
	int selectAlgorithm, int selectIntegralEqu,
	const std::string& selectMono_Dual, const std::string& polarization,
	double inc_th, double inc_ph,
	double sca_th_s, double sca_ph_s, double sca_th_f, double sca_ph_f,
	double step, int N_points) {
	RCSExportConfig cfg(inputFile);

	if (selectAlgorithm == 0) {
		cfg.selectAlgo = selectAlgorithm::MoM;
	}
	else if (selectAlgorithm == 1) {
		cfg.selectAlgo = selectAlgorithm::FMM;
	}
	else if (selectAlgorithm == 2) {
		cfg.selectAlgo = selectAlgorithm::MLFMA;
	}
	else {
		throw std::invalid_argument("Invalid selectAlgorithm (must be 0,1,2)");
	}

	cfg.material = materialType::PEC;

	cfg.integralEqu = (selectIntegralEqu == 0) ? IntegralEquation::EFIE : IntegralEquation::CFIE;

	cfg.monoDual = (selectMono_Dual == "mono") ? MonoDual::Mono : MonoDual::Dual;

	cfg.polWave = (polarization == "v") ? Polarization::V : Polarization::H;

	cfg.nPoints = N_points;

	cfg.incTheta = inc_th;
	cfg.incPhi = inc_ph;

	cfg.scaThetaStart = sca_th_s;
	cfg.scaThetaEnd = sca_th_f;
	cfg.scaPhiStart = sca_ph_s;
	cfg.scaPhiEnd = sca_ph_f;
	cfg.scaStep = step;

	// cfg.setFixedTimeStamp("...");

	cfg.validate();

	return cfg;
}

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
			BatchRunner::run(argv[2]);
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
	const std::string inputFile = "D:\\MyCode\\PMCHWT_MLFMA\\DATA\\Cube_0d5_0f20\\Cube_0d5_0f20_1e9.nas";
	const std::string movingPathFile = "D:\\MyCode\\PMCHWT_MLFMA\\DATA\\Cube_0d5_0f20\\moving_path_Cube_1e9_z.txt";
	/************************************************************************/
	double freq = 1.0e9, E0 = 1.0;
	double inc_th = 90.0, inc_ph = 0.0;
	double sca_th_s = 90.0, sca_th_f = 90.0;
	double sca_ph_s = 0.0, sca_ph_f = 360.0;
	double step = 1;
	int N_points = 7;// 1, 3, 4, 6, 7, 12, 13, 16

	std::string selectMono_Dual = "mono";// mono / dual
	std::string polarization = "v";// horizontal->90 / vertical->0
	int selectAlgorithm = 2;// 0==>MoM;1==>FMM;2==>MLFMM
	int selectIntegralEqu = 1;// 0==>EFIE;1==>CFIE;2==>PMCHWT
	int selectMatrixSolver = 0;// 0==>GMRES;1==>CGS
	// ========= Mono-FSB-MLFMA =========
	bool useFSBMono = false;// true==>open FSBMono;false==>off
	double FSBMonoFactor = 1.2;// FSB Ultra-bandwidth factor==>rho
	bool exportFSBComplexField = false;// Export the complex electric field after FSB interpolation
	bool exportFSBSamples = false;// Export the FSB samples after FSB interpolation
	// ======= Radar Position ========
	double radarX = 0.5; // monostatic radar global x
	double radarY = 0.5;
	double radarZ = 0.7071067812;
	// ======= FSB-MLFMA-Moving ========
	bool useFSBMLFMAMoving = false; // FSB moving path mode
	double movingAngleMarginDeg = 2.0; // FSB moving angular margin
	double movingAngleMergeTolDeg = 0.25; // moving angle merge tolerance
	bool exportMovingAngles = true; // output moving angle map
	bool exportMovingComplexField = false; // output moving complex field
	bool exportMovingSamples = true; // output FSB moving samples
	// ======= Accurate-MLFMA-Moving ========
	bool useAccurateMLFMAMoving = true; // true ==> solve every moving state accurately
	// ======= ACA-RHS-MLFMA-Moving ========
	bool useACARHSMoving = false; // true ==> open ACA-RHS-MLFMA-Moving
	double acaRHSTol = 1.0e-3;// ACA low-rank approximation tolerance
	int acaRHSMaxRank = 0;// ACA maximum rank constraint, 0==>No manual restriction; >0==>Use up to the specified rank
	double acaRHSPivotEps = 1.0e-14;// ACA Pivot Threshold
	int acaRHSCheckCount = 10;// ACA inspection sample size
	bool exportACARHSCheck = true;// export ACA decomposition factor
	bool exportACARHSFactors = false;// Factor matrices obtained from ACA-RHS low-rank decomposition
	// ======= Die =======
	std::complex<double> epsilonR(1.0, -0.01);
	std::complex<double> muR(1.0, 0.0);
	/************************************************************************/
	omp_set_dynamic(0);
	omp_set_num_threads(12);
	/************************************************************************/
	std::ios::sync_with_stdio(false);
	std::cin.tie(nullptr);
	/************************************************************************/
	int debugInfoToFile = 0;// 0=>off�?=>on
	/************************************************************************/
	RCSExportConfig cfg = createRCSExportConfig(inputFile, selectAlgorithm, selectIntegralEqu, selectMono_Dual, polarization,
		inc_th, inc_ph, sca_th_s, sca_ph_s, sca_th_f, sca_ph_f, step, N_points);
	cfg.useFSBMono = useFSBMono;
	cfg.fsbRho = FSBMonoFactor;
	cfg.exportFSBComplexField = exportFSBComplexField;
	cfg.exportFSBSamples = exportFSBSamples;
	cfg.useFSBMLFMAMoving = useFSBMLFMAMoving;
	cfg.movingPathFile = movingPathFile;
	cfg.radarX = radarX;
	cfg.radarY = radarY;
	cfg.radarZ = radarZ;
	cfg.movingAngleMarginDeg = movingAngleMarginDeg;
	cfg.movingAngleMergeTolDeg = movingAngleMergeTolDeg;
	cfg.exportMovingAngles = exportMovingAngles;
	cfg.exportMovingComplexField = exportMovingComplexField;
	cfg.exportMovingSamples = exportMovingSamples;
	cfg.useAccurateMLFMAMoving = useAccurateMLFMAMoving;
	cfg.useACARHSMoving = useACARHSMoving;
	cfg.acaRHSTol = acaRHSTol;
	cfg.acaRHSMaxRank = acaRHSMaxRank;
	cfg.acaRHSPivotEps = acaRHSPivotEps;
	cfg.exportACARHSCheck = exportACARHSCheck;
	cfg.acaRHSCheckCount = acaRHSCheckCount;
	cfg.exportACARHSFactors = exportACARHSFactors;
	/************************************************************************/
	double polVal = (polarization == "v") ? 0.0 : 90.0;
	EMSource wave(freq, polVal);

	if (selectIntegralEqu == 2) {
		cfg.material = materialType::Die;
		cfg.integralEqu = IntegralEquation::PMCHWT;
		cfg.epsilonR = epsilonR;
		cfg.muR = muR;
		wave.initDielectric(epsilonR, muR);
	}

	cfg.validate();
	ScopedRunLogger runLogger(cfg.buildFilePath("_RUN_LOG", "txt"));
	if (runLogger.enabled()) {
		std::cout << "Run log file: " << runLogger.path() << std::endl;
	}

	gaussPoints gausspoint(N_points);

	if (readNasData(inputFile, points, faces) != 0) {
		std::cerr << "Error: Failed to read Nastran file (main.cpp)" << std::endl;
		return 1;
	}

	std::vector<RWGBase> rwgBases;
	if (!RWGGenerator::GenerateRWG(faces, points, rwgBases)) {
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

	auto startOCTree = std::chrono::high_resolution_clock::now();
	OCTree octree(rwgBases, wave.wavelength(), inputFile);
	auto endOCTree = std::chrono::high_resolution_clock::now();
	auto durationOCTree = std::chrono::duration_cast<std::chrono::microseconds>(endOCTree - startOCTree);
	std::cout << "<<<< Info: OCTree built, elapsed: " << durationOCTree.count() / 1000.0 << " ms >>>> (main.cpp)\n" << std::endl;

	if (selectMatrixSolver != 0 && selectMatrixSolver != 1) {
		std::cerr << "Error: Invalid selectMatrixSolver value: "
			<< selectMatrixSolver
			<< " (0 = GMRES, 1 = CGS)" << std::endl;
		return 3;
	}

	if (selectAlgorithm == 0) {
		auto startMoM = std::chrono::high_resolution_clock::now();
		MoM mom(cfg, selectIntegralEqu, selectMatrixSolver, selectMono_Dual, polarization,
			octree.octreeNodes_, rwgBases, octree.maxLevel_, gausspoint,
			E0, wave);
		auto endMoM = std::chrono::high_resolution_clock::now();
		auto durationMoM = std::chrono::duration_cast<std::chrono::microseconds>(endMoM - startMoM);
		std::cout << "<<<< Info: MoM done, elapsed: " << durationMoM.count() / 1000.0 << " ms >>>> (main.cpp)\n" << std::endl;
	}
	else if (selectAlgorithm == 1) {
		auto startFMM = std::chrono::high_resolution_clock::now();
		FMM fmm(cfg, selectIntegralEqu, selectMatrixSolver, selectMono_Dual, polarization,
			octree.octreeNodes_, rwgBases, octree.maxLevel_, gausspoint,
			E0, wave);
		auto endFMM = std::chrono::high_resolution_clock::now();
		auto durationFMM = std::chrono::duration_cast<std::chrono::microseconds>(endFMM - startFMM);
		std::cout << "<<<< Info: FMM done, elapsed: " << durationFMM.count() / 1000.0 << " ms >>>> (main.cpp)\n" << std::endl;
	}
	else if (selectAlgorithm == 2) {
		auto startMLFMA = std::chrono::high_resolution_clock::now();
		MLFMM mlfmm(cfg, selectIntegralEqu, selectMatrixSolver, selectMono_Dual, polarization,
			octree.octreeNodes_, octree.octreeNodesDRvec_, rwgBases, octree.maxLevel_, gausspoint,
			E0, wave);
		auto endMLFMA = std::chrono::high_resolution_clock::now();
		auto durationMLFMA = std::chrono::duration_cast<std::chrono::microseconds>(endMLFMA - startMLFMA);
		std::cout << "<<<< Info: MLFMA done, elapsed: " << durationMLFMA.count() / 1000.0 << " ms >>>> (main.cpp)\n" << std::endl;
	}
	else {
		std::cerr << "Error: Invalid selectAlgorithm value: " << selectAlgorithm << std::endl;
	}

	return 0;
}