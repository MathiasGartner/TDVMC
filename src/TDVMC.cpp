#include "mpi.h"

#include "ConfigItem.h"
#include "Constants.h"
#include "ICorrelatedSamplingData.h"
#include "MathOperators.h"
#include "MPIMethods.h"
#include "SimulationStepData.h"
#include "Timer.h"
#include "Utils.h"
#include "VMCSampler.h"

#include "PhysicalSystems/IPhysicalSystem.h"
#include "PhysicalSystems/BosonsBulk.h"
#include "PhysicalSystems/BosonsBulkDamped.h"
#include "PhysicalSystems/BosonCluster.h"
#include "PhysicalSystems/BosonClusterWithLog.h"
#include "PhysicalSystems/BosonClusterWithLogParam.h"
#include "PhysicalSystems/BosonMixtureCluster.h"
#include "PhysicalSystems/BulkSplines.h"
#include "PhysicalSystems/BulkSplinesPhi.h"
#include "PhysicalSystems/BulkSplinesScaled.h"
#include "PhysicalSystems/BulkQT.h"
#include "PhysicalSystems/BulkQTPhi.h"
#include "PhysicalSystems/GaussianWavepacket.h"
#include "PhysicalSystems/HardSphereBosons.h"
#include "PhysicalSystems/HardSphereBosonsExp.h"
#include "PhysicalSystems/HeBulk.h"
#include "PhysicalSystems/HeDrop.h"
#include "PhysicalSystems/NUBosonsBulk.h"
#include "PhysicalSystems/PBosonsBulk.h"

#include "Potentials/PotentialManager.h"

#include "test/MPITest.h"
#include "test/PiCalculator.h"
#include "test/Tests.h"

//#undef SEEK_SET
//#undef SEEK_END
//#undef SEEK_CUR

//#include <armadillo>
#include <cstring>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <math.h>
#include <signal.h>
#include <sstream>
#include <stdlib.h>
#include <string>
#include <unistd.h>
#include <vector>

//#include "../resources/armadillo"
#include "../resources/json/json-forwards.h"
#include "../resources/json/json.h"

using namespace std;

PhysicalSystems::IPhysicalSystem* sys;

vector<ConfigItem> configItems;

string CONFIG_VERSION;
string SYSTEM_TYPE;
string OUTPUT_DIRECTORY;		//from config file
string configDirectory;
string configFilePath;
string OUT_DIR;					//directory name generated from parameter settings
int N;           	    		//number of particles
int DIM;     	        	 	//number of dimensions
int N_PARAM;	        	  	//number of parameters of trial function
double MC_STEP;
double MC_STEP_OFFSET;			//INFO: currently not used
int MC_NSTEPS;
int MC_NTHERMSTEPS;
int MC_NINITIALIZATIONSTEPS;
int MC_VERY_FIRST_NINITIALIZATIONSTEPS;
int MC_NADDITIONALSTEPS;
int MC_NADDITIONALTHERMSTEPS;
int MC_NADDITIONALINITIALIZATIONSTEPS;
double RHO;
double RC;          			//cutoff for WF and LJ
double TIMESTEP;
double TOTALTIME;
int IMAGINARY_TIME;
int ODE_SOLVER_TYPE;
int USE_PRECONDITIONING;
double LBOX;
double LBOX_R;
vector<double> PARAMS_REAL;
vector<double> PARAMS_IMAGINARY;
double PARAM_PHIR;
double PARAM_PHII;
int USE_PARAMETER_ACCEPTANCE_CHECK;
int PARAMETER_ACCEPTANCE_CHECK_TYPE;
int WRITE_EVERY_NTH_STEP_TO_FILE;
int WRITE_SINGLE_FILES;
int MC_NSTEP_MULTIPLICATION_FACTOR_FOR_WRITE_DATA;
int USE_MEAN_FOR_FINAL_PARAMETERS;
int USE_NORMALIZE_WF;
int USE_ADJUST_PARAMETERS;
int UPDATE_SAMPLES_EVERY_NTH_STEP;
double UPDATE_SAMPLES_PERCENT;
int GR_BIN_COUNT;
int USE_NURBS;
vector<double> NURBS_GRID;
vector<int> PARTICLE_TYPES;
vector<double> SYSTEM_PARAMS;

string requiredConfigVersion = "0.18";
int numOfProcesses = 1;
int rootRank = 0;
int processRank = 0;
long long nAcceptances = 0;
long long nTrials = 0;
double mc_nsteps;
int mc_nsteps_original;
double mc_nadditionalsteps;
double currentTime;

vector<double> localOperators; // for <O_k>
double localEnergyR; // for <E^R>
double localEnergyI; // for <E^I>
vector<vector<double> > localOperatorsMatrix; // for <O_k O_j>
vector<double> localOperatorlocalEnergyR; // for <O_k E^R>
vector<double> localOperatorlocalEnergyI; // for <O_k E^I>
vector<double> otherExpectationValues; // eg. for potential and kinetic energy
vector<double> additionalSystemProperties; // for properties at the end of the simulation

bool doNotAcceptStep = false;
vector<vector<double> > uRList;
vector<vector<double> > uIList;
vector<double> phiRList;
vector<double> phiIList;
vector<vector<double> > uRListDiffs;
vector<vector<double> > uIListDiffs;
vector<double> phiRListDiffs;
vector<double> phiIListDiffs;

vector<double> AllLocalEnergyR;
vector<double> AllLocalEnergyI;
vector<vector<double> > AllLocalOperators;
vector<vector<double> > AllOtherExpectationValues;
vector<vector<double> > AllParametersR;
vector<vector<double> > AllParametersI;
vector<vector<double> > AllAdditionalSystemProperties;

int currentSampleIndexForUpdate;
vector<ICorrelatedSamplingData*> correlatedSamplingData;

vector<double> times;
vector<double> previousStepWeights = { 0.982014, 0.952574, 0.880797, 0.731059, 0.5, 0.268941, 0.119203, 0.0474259, 0.0179862, 0.00669285 };
vector<SimulationStepData> previousStepData;

mt19937_64 generator;
//default_random_engine generator;
uniform_real_distribution<double> distUniform(0.0, 1.0);
uniform_int_distribution<int> distParticleIndex;
normal_distribution<double> distNormal(0.0, 1.0);

////////////////////
/// Log messages ///
////////////////////

enum MessageType
{
	NORMAL, WARNING, ERROR
};

void Log(string message)
{
	cout << "#" << setfill(' ') << setw(2) << processRank << "@" << get_cpu_id() << ": " << message << endl << flush;
}

void Log(string message, MessageType messageType)
{
	switch (messageType)
	{
	case NORMAL:
		Log(message);
		break;
	case WARNING:
		cout << "\033[1;33m" << "#" << setfill(' ') << setw(2) << processRank << "@" << get_cpu_id() << ": " << message << "\033[0m" << endl << flush;
		break;
	case ERROR:
		cout << "\033[1;31m" << "#" << setfill(' ') << setw(2) << processRank << "@" << get_cpu_id() << ": " << message << "\033[0m" << endl << flush;
		break;
	}
}

void PrintData(vector<double>& data)
{
	for (unsigned int i = 0; i < data.size(); i++)
	{
		cout << data[i] << ", ";
	}
	cout << endl;
}

void PrintData(vector<bool>& data)
{
	for (unsigned int i = 0; i < data.size(); i++)
	{
		cout << (data[i] ? '#' : '_');
		if (i % 100 == 99)
		{
			cout << endl;
		}
	}
	cout << endl;
}

void PrintData(vector<vector<double> > data)
{
	for (unsigned int i = 0; i < data.size(); i++)
	{
		PrintData(data[i]);
	}
	cout << endl;
}

///////////////////////
/// signal handling ///
///////////////////////

void onSignalStop(int signum)
{
	if (signum == 10) //SIGUSR1
	{
		cout << "Received SIGUSR1!" << endl;
		cout << "Finishing simulation at t=" << currentTime << endl;
		TOTALTIME = currentTime;
	}
}

////////////////////////////////
/// random number generation ///
////////////////////////////////

double random01()
{
	return distUniform(generator);
}

double randomNormal()
{
	return distNormal(generator);
}

double randomNormal(double sigma)
{
	return randomNormal() * sigma;
}

double randomNormal(double sigma, double mu)
{
	return randomNormal() * sigma + mu;
}

int randomInt(int maxValue)
{
	//INFO: test for maxValue > 0 is only needed for Gaussian wavepacket simulation where only one  particle is used
	return maxValue > 0 ? ((int) floor(random01() * maxValue)) % maxValue : maxValue;
	//return ((int) floor(random01() * maxValue)) % maxValue;
}

int randomInt(int minValue, int maxValue)
{
	return randomInt(maxValue - minValue) + minValue;
}

int randomParticleIndex()
{
	return distParticleIndex(generator);
}

void ReadRandomGeneratorStatesFromFile(string fileNamePrefix)
{
	ifstream fGenerator(configDirectory + fileNamePrefix + "_generator_" + to_string(processRank) + ".dat");
	fGenerator >> generator;
	fGenerator.close();
	ifstream fUniform(configDirectory + fileNamePrefix + "_uniform_" + to_string(processRank) + ".dat");
	fUniform >> distUniform;
	fUniform.close();
	ifstream fNormal(configDirectory + fileNamePrefix + "_normal_" + to_string(processRank) + ".dat");
	fNormal >> distNormal;
	fNormal.close();
	ifstream fParticle(configDirectory + fileNamePrefix + "_particleIndex_" + to_string(processRank) + ".dat");
	fParticle >> distParticleIndex;
	fParticle.close();
}

void WriteRandomGeneratorStatesToFile(string fileNamePrefix)
{
	ofstream fGenerator(OUT_DIR + fileNamePrefix + "_generator_" + to_string(processRank) + ".dat");
	fGenerator << generator;
	fGenerator.close();
	ofstream fUniform(OUT_DIR + fileNamePrefix + "_uniform_" + to_string(processRank) + ".dat");
	fUniform << distUniform;
	fUniform.close();
	ofstream fNormal(OUT_DIR + fileNamePrefix + "_normal_" + to_string(processRank) + ".dat");
	fNormal << distNormal;
	fNormal.close();
	ofstream fParticle(OUT_DIR + fileNamePrefix + "_particleIndex_" + to_string(processRank) + ".dat");
	fParticle << distParticleIndex;
	fParticle.close();
}

void (*RandomDisplaceParticle_DIM)(vector<double>& r);

void RandomDisplaceParticle_1D(vector<double>& r)
{
	r[0] += randomNormal(MC_STEP);
}

void RandomDisplaceParticle_2D(vector<double>& r)
{
	r[0] += randomNormal(MC_STEP);
	r[1] += randomNormal(MC_STEP);
}

void RandomDisplaceParticle_3D(vector<double>& r)
{
	r[0] += randomNormal(MC_STEP);
	r[1] += randomNormal(MC_STEP);
	r[2] += randomNormal(MC_STEP);
}

void RandomDisplaceParticle(vector<double>& r)
{
	for (unsigned int i = 0; i < r.size(); i++)
	{
		r[i] += randomNormal(MC_STEP);
	}
}

///////////////////////////
/// configuration items ///
///////////////////////////

void RegisterAllConfigItems()
{
	configItems.push_back(ConfigItem("CONFIG_VERSION", &CONFIG_VERSION, ConfigItemType::STRING));
	configItems.push_back(ConfigItem("SYSTEM_TYPE", &SYSTEM_TYPE, ConfigItemType::STRING));
	configItems.push_back(ConfigItem("OUTPUT_DIRECTORY", &OUTPUT_DIRECTORY, ConfigItemType::STRING));
	configItems.push_back(ConfigItem("N", &N, ConfigItemType::INT));
	configItems.push_back(ConfigItem("DIM", &DIM, ConfigItemType::INT));
	configItems.push_back(ConfigItem("N_PARAM", &N_PARAM, ConfigItemType::INT));
	configItems.push_back(ConfigItem("RHO", &RHO, ConfigItemType::DOUBLE));
	configItems.push_back(ConfigItem("RC", &RC, ConfigItemType::DOUBLE));
	configItems.push_back(ConfigItem("MC_STEP", &MC_STEP, ConfigItemType::DOUBLE));
	configItems.push_back(ConfigItem("MC_STEP_OFFSET", &MC_STEP_OFFSET, ConfigItemType::DOUBLE));
	configItems.push_back(ConfigItem("MC_NSTEPS", &MC_NSTEPS, ConfigItemType::INT));
	configItems.push_back(ConfigItem("MC_NTHERMSTEPS", &MC_NTHERMSTEPS, ConfigItemType::INT));
	configItems.push_back(ConfigItem("MC_NINITIALIZATIONSTEPS", &MC_NINITIALIZATIONSTEPS, ConfigItemType::INT));
	configItems.push_back(ConfigItem("MC_VERY_FIRST_NINITIALIZATIONSTEPS", &MC_VERY_FIRST_NINITIALIZATIONSTEPS, ConfigItemType::INT));
	configItems.push_back(ConfigItem("MC_NADDITIONALSTEPS", &MC_NADDITIONALSTEPS, ConfigItemType::INT));
	configItems.push_back(ConfigItem("MC_NADDITIONALTHERMSTEPS", &MC_NADDITIONALTHERMSTEPS, ConfigItemType::INT));
	configItems.push_back(ConfigItem("MC_NADDITIONALINITIALIZATIONSTEPS", &MC_NADDITIONALINITIALIZATIONSTEPS, ConfigItemType::INT));
	configItems.push_back(ConfigItem("TIMESTEP", &TIMESTEP, ConfigItemType::DOUBLE));
	configItems.push_back(ConfigItem("TOTALTIME", &TOTALTIME, ConfigItemType::DOUBLE, true));
	configItems.push_back(ConfigItem("IMAGINARY_TIME", &IMAGINARY_TIME, ConfigItemType::INT));
	configItems.push_back(ConfigItem("ODE_SOLVER_TYPE", &ODE_SOLVER_TYPE, ConfigItemType::INT, true));
	configItems.push_back(ConfigItem("USE_PRECONDITIONING", &USE_PRECONDITIONING, ConfigItemType::INT, true));
	configItems.push_back(ConfigItem("USE_PARAMETER_ACCEPTANCE_CHECK", &USE_PARAMETER_ACCEPTANCE_CHECK, ConfigItemType::INT, true));
	configItems.push_back(ConfigItem("PARAMETER_ACCEPTANCE_CHECK_TYPE", &PARAMETER_ACCEPTANCE_CHECK_TYPE, ConfigItemType::INT, true));
	configItems.push_back(ConfigItem("WRITE_EVERY_NTH_STEP_TO_FILE", &WRITE_EVERY_NTH_STEP_TO_FILE, ConfigItemType::INT, true));
	configItems.push_back(ConfigItem("WRITE_SINGLE_FILES", &WRITE_SINGLE_FILES, ConfigItemType::INT, true));
	configItems.push_back(ConfigItem("MC_NSTEP_MULTIPLICATION_FACTOR_FOR_WRITE_DATA", &MC_NSTEP_MULTIPLICATION_FACTOR_FOR_WRITE_DATA, ConfigItemType::INT, true));
	configItems.push_back(ConfigItem("USE_MEAN_FOR_FINAL_PARAMETERS", &USE_MEAN_FOR_FINAL_PARAMETERS, ConfigItemType::INT, true));
	configItems.push_back(ConfigItem("USE_NORMALIZE_WF", &USE_NORMALIZE_WF, ConfigItemType::INT, true));
	configItems.push_back(ConfigItem("USE_ADJUST_PARAMETERS", &USE_ADJUST_PARAMETERS, ConfigItemType::INT, true));
	configItems.push_back(ConfigItem("UPDATE_SAMPLES_EVERY_NTH_STEP", &UPDATE_SAMPLES_EVERY_NTH_STEP, ConfigItemType::INT, true));
	configItems.push_back(ConfigItem("UPDATE_SAMPLES_PERCENT", &UPDATE_SAMPLES_PERCENT, ConfigItemType::DOUBLE, true));
	configItems.push_back(ConfigItem("GR_BIN_COUNT", &GR_BIN_COUNT, ConfigItemType::INT));
	configItems.push_back(ConfigItem("USE_NURBS", &USE_NURBS, ConfigItemType::INT));
	configItems.push_back(ConfigItem("NURBS_GRID", NURBS_GRID, ConfigItemType::ARR_DOUBLE));
	configItems.push_back(ConfigItem("PARTICLE_TYPES", PARTICLE_TYPES, ConfigItemType::ARR_INT));
	configItems.push_back(ConfigItem("SYSTEM_PARAMS", SYSTEM_PARAMS, ConfigItemType::ARR_DOUBLE));
	configItems.push_back(ConfigItem("PARAMS_REAL", PARAMS_REAL, ConfigItemType::ARR_DOUBLE));
	configItems.push_back(ConfigItem("PARAMS_IMAGINARY", PARAMS_IMAGINARY, ConfigItemType::ARR_DOUBLE));
	configItems.push_back(ConfigItem("PARAM_PHIR", &PARAM_PHIR, ConfigItemType::DOUBLE));
	configItems.push_back(ConfigItem("PARAM_PHII", &PARAM_PHII, ConfigItemType::DOUBLE));
}

void ReadConfig(string filePath)
{
	//cout << "file exists: " << FileExist(filePath) << "." << endl;

	string errors;
	bool successful;
	Json::Value configData;
	Json::CharReaderBuilder builder;
	ifstream configFile(filePath, ifstream::binary);

	successful = Json::parseFromStream(builder, configFile, &configData, &errors);

	if (successful)
	{
		for (auto ci : configItems)
		{
			//cout << "Read: " << ci.name << endl;
			ci.setValue(configData[ci.name]);
		}
	}
	else
	{
		Log("could not read config file", ERROR);
	}
}

void PrintConfig()
{
	Log("===== Configuration ====");
	for (auto ci : configItems)
	{
		cout << ci.name << ": " << ci.getJsonString() << endl;
	}
	cout << "=============================" << endl << endl;
}

void WriteConfig(string fileName, vector<double>& uR, vector<double>& uI, double phiR, double phiI)
{
	ofstream configFile;
	configFile.open(OUT_DIR + fileName, ios::out);
	configFile.precision(8);

	configFile << "{" << endl;
	for (unsigned int i = 0; i < configItems.size(); i++)
	{
		auto ci = configItems[i];
		configFile << "\t" << "\"" << ci.name << "\"" << " : " << ci.getJsonString();
		if (i == configItems.size() - 1)
		{
			configFile << endl;
		}
		else
		{
			configFile << "," << endl;
		}
	}
	configFile << "}";

	configFile.close();
}

void BroadcastConfigItem(ConfigItem& ci)
{
	if (ci.type == STRING)
	{
		MPIMethods::BroadcastValue(ci.variableString, 300);
	}
	else if (ci.type == INT)
	{
		MPIMethods::BroadcastValue(ci.variableInt);
	}
	else if (ci.type == DOUBLE)
	{
		MPIMethods::BroadcastValue(ci.variableDouble);
	}
	else if (ci.type == ARR_INT)
	{
		int size = ci.variableArrInt->size();
		MPIMethods::BroadcastValue(&size);
		if (processRank != rootRank)
		{
			ci.variableArrInt->resize(size);
		}
		MPIMethods::BroadcastValues(*(ci.variableArrInt));
	}
	else if (ci.type == ARR_DOUBLE)
	{
		int size = ci.variableArrDouble->size();
		MPIMethods::BroadcastValue(&size);
		if (processRank != rootRank)
		{
			ci.variableArrDouble->resize(size);
		}
		MPIMethods::BroadcastValues(*(ci.variableArrDouble));
	}
}

void BroadcastChangeableConfigItems()
{
	for (auto ci : configItems)
	{
		if (ci.allowChangeAtRuntime)
		{
			BroadcastConfigItem(ci);
		}
	}
}

void BroadcastConfig()
{
	for (auto ci : configItems)
	{
		//cout << "Broadcast: " << ci.name << endl;
		BroadcastConfigItem(ci);
	}
}

void CreateOutputDirectory()
{
	//INFO: append config options to output directory path, create the directory and copy the config file
	int tmp;
	ostringstream strs;
	strs << "step=" << MC_NSTEPS << "_therm=" << MC_NTHERMSTEPS << "_time=" << TIMESTEP;
	OUT_DIR = OUTPUT_DIRECTORY + strs.str() + "/";
	tmp = system(("rm -rf " + OUT_DIR).c_str());
	tmp = system(("mkdir " + OUT_DIR).c_str());
	CopyFile(configFilePath, OUT_DIR + "vmc.config");
	cout << tmp << endl;
}

void WriteParticleInputFile(string fileName, vector<vector<double> >& R)
{
	vector<vector<double> > particlePositionList = { { } };
	for (int i = 0; i < N; i++)
	{
		for (int a = 0; a < DIM; a++)
		{
			particlePositionList[0].push_back(R[i][a]);
		}
	}
	WriteDataToFile(particlePositionList, fileName, "", 1, false);
}

void BroadcastNewParameters(vector<double>& uR, vector<double>& uI, double* phiR, double* phiI)
{
	MPIMethods::BroadcastValues(uR);
	MPIMethods::BroadcastValues(uI);
	MPIMethods::BroadcastValue(phiR);
	MPIMethods::BroadcastValue(phiI);
}

void Init()
{
	if (FileExist(configDirectory + "random/state" + "_generator_" + to_string(processRank) + ".dat") && FileExist(configDirectory + "random/state" + "_uniform_" + to_string(processRank) + ".dat") && FileExist(configDirectory + "random/state" + "_normal_" + to_string(processRank) + ".dat"))
	{
		Log("read random generator configurations");
		ReadRandomGeneratorStatesFromFile("random/state");
	}
	else
	{
		generator = mt19937_64(processRank + 1);
		//generator = default_random_engine(processRank + 1);
		distParticleIndex = uniform_int_distribution<int>(0, N - 1);
		srand(processRank + 1);
	}
	LBOX = pow((N / RHO), 1.0 / DIM); //box with dimensions [-L/2, L/2]
	LBOX_R = 1.0 / LBOX;
	nAcceptances = 0;
	nTrials = 0;

	mc_nsteps = (double) MC_NSTEPS;
	mc_nsteps_original = MC_NSTEPS;
	mc_nadditionalsteps = (double) MC_NADDITIONALSTEPS;

	currentSampleIndexForUpdate = 0;
	correlatedSamplingData.resize(MC_NSTEPS);
	//TODO: Let the IPhysicalSystem sys create the objects needed for storing the correlated sampling data
	for (int i = 0; i < MC_NSTEPS; i++)
	{
		correlatedSamplingData[i] = new CSDataBulkSplines();
	}

	//INFO: Init Utils and dimension-dependent functions
	switch (DIM)
	{
	case 1:
		VectorDotProduct_DIM = VectorDotProduct_1D;
		VectorNorm2_DIM = VectorNorm2_1D;
		VectorNorm_DIM = VectorNorm_1D;
		VectorDisplacement_DIM = VectorDisplacement_1D;
		GetVectorNIC_DIM = GetVectorNIC_1D;
		VectorDiffNIC_DIM = VectorDiffNIC_1D;
		VectorDisplacementNIC_DIM = VectorDisplacementNIC_1D;
		AssignVector_DIM = AssignVector_1D;
		RandomDisplaceParticle_DIM = RandomDisplaceParticle_1D;
		break;
	case 2:
		VectorDotProduct_DIM = VectorDotProduct_2D;
		VectorNorm2_DIM = VectorNorm2_2D;
		VectorNorm_DIM = VectorNorm_2D;
		VectorDisplacement_DIM = VectorDisplacement_2D;
		GetVectorNIC_DIM = GetVectorNIC_2D;
		VectorDiffNIC_DIM = VectorDiffNIC_2D;
		VectorDisplacementNIC_DIM = VectorDisplacementNIC_2D;
		AssignVector_DIM = AssignVector_2D;
		RandomDisplaceParticle_DIM = RandomDisplaceParticle_2D;
		break;
	case 3:
		VectorDotProduct_DIM = VectorDotProduct_3D;
		VectorNorm2_DIM = VectorNorm2_3D;
		VectorNorm_DIM = VectorNorm_3D;
		VectorDisplacement_DIM = VectorDisplacement_3D;
		GetVectorNIC_DIM = GetVectorNIC_3D;
		VectorDiffNIC_DIM = VectorDiffNIC_3D;
		VectorDisplacementNIC_DIM = VectorDisplacementNIC_3D;
		AssignVector_DIM = AssignVector_3D;
		RandomDisplaceParticle_DIM = RandomDisplaceParticle_3D;
		break;
	default:
		VectorDotProduct_DIM = VectorDotProduct;
		VectorNorm2_DIM = VectorNorm2;
		VectorNorm_DIM = VectorNorm;
		VectorDisplacement_DIM = VectorDisplacement;
		GetVectorNIC_DIM = GetVectorNIC;
		VectorDiffNIC_DIM = VectorDiffNIC;
		VectorDisplacementNIC_DIM = VectorDisplacementNIC;
		AssignVector_DIM = AssignVector;
		RandomDisplaceParticle_DIM = RandomDisplaceParticle;
		break;
	}
}

void PostSystemInit()
{
	localOperators.resize(N_PARAM);
	localOperatorsMatrix.resize(N_PARAM);
	for (auto &row : localOperatorsMatrix)
	{
		row.resize(N_PARAM);
	}
	localOperatorlocalEnergyR.resize(N_PARAM);
	localOperatorlocalEnergyI.resize(N_PARAM);
	otherExpectationValues.resize(sys->GetNumOfOtherExpectationValues());
	additionalSystemProperties.resize(sys->GetNumOfAdditionalSystemProperties());
	AllAdditionalSystemProperties.resize(0);
}

void WriteParticlesToFile(vector<vector<double> >& R, string ending)
{
	WriteDataToFile(R, "position" + ending, "x, y, z");
}

bool LoadLastPositionsFromFile(string filename, vector<vector<double> >& R)
{
	bool successful = false;
	string line;
	string prevline;
	vector<string> coordinates;
	ifstream file;
	//Log("try read particle configuration from file: " + configDirectory + filename + ".csv");
	file.open(configDirectory + filename + ".csv", ios::in);
	while (getline(file, line))
	{
		prevline = line;
	}
	if (prevline.length() > 0)
	{
		Log("init coordinates from file: " + filename);
		coordinates = split(prevline, ',');
		int j = 0;
		for (unsigned int i = 0; i < coordinates.size(); i += DIM)
		{
			for (int a = 0; a < DIM; a++)
			{
				R[j][a] = stod(coordinates[i + a]);
			}
			j++;
		}
		successful = true;
	}
	file.close();
	return successful;
}

void InitCoordinateConfiguration(vector<vector<double> >& R)
{
	int fileFound = 0;
	if (processRank == rootRank)
	{
		string filename = "coords/particleconfiguration_" + to_string(N) + "_" + to_string(DIM) + "D_" + to_string(processRank);
		if (LoadLastPositionsFromFile(filename, R))
		{
			fileFound = 1;
		}
	}
	MPIMethods::BroadcastValue(&fileFound);
	if (fileFound)
	{
		//MPIMethods::BroadcastValues(R);
		string filename = "coords/particleconfiguration_" + to_string(N) + "_" + to_string(DIM) + "D_" + to_string(processRank);
		LoadLastPositionsFromFile(filename, R);
	}
	else
	{
		if (processRank == rootRank)
		{
			cout << "LBOX=" << LBOX << endl;
		}
		int type = 2; //INFO: 1: lattice, 2: drop
		if (type == 0)
		{
			//Random
			Log("init coordinates randomly");
			for (int i = 0; i < N; i++)
			{
				for (int j = 0; j < DIM; j++)
				{
					R[i][j] = (random01() - 0.5) * LBOX / 2.0;
					if (j == 1 || j == 2)
					{
						//R[i][j] = 0;
					}
				}
			}
		}
		else if (type == 1)
		{
			//Lattice
			//TODO: funktioniert nicht
			Log("init coordinates on lattice");
			int n = round(pow(N, 1.0 / ((double) DIM)));
			double l = LBOX / n;
			l *= 0.5;
			int p = 0;
			if (DIM == 3)
			{
				for (int i = 0; i < n; i++)
				{
					for (int j = 0; j < n; j++)
					{
						for (int k = 0; k < n; k++)
						{
							R[p][0] = i * l + (random01() - 0.5) * l / 10.0;
							R[p][1] = j * l + (random01() - 0.5) * l / 10.0;
							R[p][2] = k * l + (random01() - 0.5) * l / 10.0;
							p++;
						}
					}
				}
			}
			else if (DIM == 2)
			{
				for (int i = 0; i < n; i++)
				{
					for (int j = 0; j < n; j++)
					{
						R[p][0] = i * l + (random01() - 0.5) * l / 10.0;
						R[p][1] = j * l + (random01() - 0.5) * l / 10.0;
						p++;
					}
				}
			}
			else if (DIM == 1)
			{
				for (int i = 0; i < n; i++)
				{
					R[p][0] = i * l + (random01() - 0.5) * l / 10.0;
					p++;
				}
			}
		}
		else if (type == 2)
		{
			//Drop
			Log("init coordinates for drop");
			int i = 0;
			double l = pow(LBOX, 1.0 / 3.0);
			vector<int> positions;
			while (i < N)
			{
				positions =
				{	0, 0, 0};
				for (int c = 0; c < i; c++)
				{
					positions[(c % 3 + i % 3) % 3]++;
				}
				for (int a = 0; a < DIM; a++)
				{
					R[i][a] = positions[a] * l + (random01() - 0.5) * l / 2.0;
				}
				i++;
			}
		}
	}
}

void MoveCoordinatesToFirstCell(vector<vector<double> >& R)
{
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < DIM; j++)
		{
			R[i][j] = GetCoordinateNIC(R[i][j]);
		}
	}
}

void MoveCenterOfMassToZero(vector<vector<double> >& R)
{
	vector<double> com;
	com = sys->GetCenterOfMass(R);
	for (int i = 0; i < N; i++)
	{
		for (int a = 0; a < DIM; a++)
		{
			R[i][a] -= com[a];
		}
	}
}

void DoMetropolisStep(vector<vector<double> >& R, vector<double>& uR, vector<double>& uI, double phiR, double phiI)
{
	//INFO: no performance increase for dimension-dependent functions measureable on asterix
	//int randomParticle = randomInt(N - 1);
	//RandomDisplaceParticle_DIM(R[randomParticle]);
	//R[randomParticle][0] += randomNormal(MC_STEP);
	//R[randomParticle][1] += randomNormal(MC_STEP);
	//AssignVector_DIM(oldPosition, R[randomParticle]);
	//R[randomParticle][0] = oldPosition[0];
	//R[randomParticle][1] = oldPosition[1];
	double p;
	bool sampleOkay = true;
	int randomParticle = randomParticleIndex();
	vector<double> oldPosition(R[randomParticle]);
	for (int i = 0; i < DIM; i++)
	{
		R[randomParticle][i] += randomNormal(MC_STEP);
	}
	double wfQuotient = sys->CalculateWFQuotient(R, uR, uI, phiR, phiI, randomParticle, oldPosition);
	//cout << wfQuotient << endl;

	//if (wfQuotient == 0 || sys->GetWfNew() == 0 || sys->GetWf() == 0)
	//{
	//	cout << "something == 0 !!!!!!!!!!!!!!!!!!!!!!!!" << endl;
	//	cout << "wfQuotient=" << wfQuotient << ", wfNew=" << Sys::wfNew << ", wf=" << Sys::wf << ", nTrials=" << nTrials << endl;
	//	sleep(1);
	//}
	if (!isfinite(wfQuotient) || !isfinite(sys->GetWfNew()) || !isfinite(sys->GetWf()))
	{
		//cout << "something != finite !!!!!!!!!!!!!!!!!!!!!!!!" << endl;
		//cout << "wfQuotient=" << wfQuotient << ", wfNew=" << sys->GetWfNew() << ", wf=" << sys->GetWf() << ", nTrials=" << nTrials << endl;
		//sleep(1);
		sampleOkay = false;
		if (!isfinite(wfQuotient) && sys->GetWfNew() > 0 && sys->GetWf() == 0) //INFO: this can happen when a new timestep is startet and the wavefunction is zero for the first particle configuration
		{
			sampleOkay = true;
			wfQuotient = 1; //INFO: accept by 100%
		}
	}

	p = random01();
	if (!sampleOkay || wfQuotient < p)
	{
		for (int i = 0; i < DIM; i++)
		{
			R[randomParticle][i] = oldPosition[i];
		}
	}
	else
	{
		sys->AcceptMove();
		nAcceptances++;
	}
	nTrials++;
}

//////////////////////////////////////////
/// strategies for correlated sampling ///
//////////////////////////////////////////

bool NeedToUpdateSamples(vector<ICorrelatedSamplingData*>& samples, vector<double>& uR, vector<double>& uI, double phiR, double phiI)
{
	bool updateNeeded = true;
	//double weight = 1.0;
	//double weightDiff = 0.0;
	//int count = samples.size();
	//int randomSample = randomInt(count - 1);
	//sys->CalculateWavefunction(samples[randomSample], uR, uI, phiR, phiI);
	////double newWf2 = sys->GetWf() * sys->GetWf();
	////weight = newWf2 / wfValues2[randomSample];
	//weight = exp(2.0 * sys->GetExponent() - 2.0 * correlatedSamplingData[randomSample].exponent);
	//weightDiff = abs(1.0 - weight);
	//updateNeeded = weightDiff > 0.99;
	return updateNeeded;
}

void UpdateSample(int sampleIndex, vector<ICorrelatedSamplingData*>& samples, vector<double>& uR, vector<double>& uI, double phiR, double phiI)
{
	sys->CalculateWavefunction(samples[sampleIndex]->R, uR, uI, phiR, phiI);
	for (int n = 0; n < MC_NTHERMSTEPS; n++)
	{
		DoMetropolisStep(samples[sampleIndex]->R, uR, uI, phiR, phiI);
	}
	sys->CalculateWavefunction(samples[sampleIndex]->R, uR, uI, phiR, phiI); //TODO: it would be sufficient to call sys->RefreshLocalOperators() since the splineSums should be kept up to date for samples[randomSample]->R during the MetropolisSteps
	sys->CalculateOtherLocalOperators(samples[sampleIndex]->R);

	samples[sampleIndex]->wf = sys->GetWf();
	samples[sampleIndex]->wf2 = samples[sampleIndex]->wf * samples[sampleIndex]->wf;
	samples[sampleIndex]->exponent = sys->GetExponent();
	samples[sampleIndex]->exponent2 = samples[sampleIndex]->exponent * samples[sampleIndex]->exponent;
	samples[sampleIndex]->localOperators = sys->GetLocalOperators();
	sys->FillCorrelatedSamplingData(samples[sampleIndex]);
}

void UpdateSamplesRandom(int nrOfSamplesToUpdate, vector<ICorrelatedSamplingData*>& samples, vector<double>& uR, vector<double>& uI, double phiR, double phiI)
{
	int count = samples.size();

	for (int i = 0; i < nrOfSamplesToUpdate; i++)
	{
		int randomSample = randomInt(count - 1);
		UpdateSample(randomSample, samples, uR, uI, phiR, phiI);
	}
}

void UpdateSamplesConsecutive(int nrOfSamplesToUpdate, vector<ICorrelatedSamplingData*>& samples, vector<double>& uR, vector<double>& uI, double phiR, double phiI)
{
	int count = samples.size();
	for (int i = 0; i < nrOfSamplesToUpdate; i++)
	{
		UpdateSample(currentSampleIndexForUpdate, samples, uR, uI, phiR, phiI);
		currentSampleIndexForUpdate = (currentSampleIndexForUpdate + 1) % count;
	}
}

void UpdateAllSamples(vector<ICorrelatedSamplingData*>& samples, vector<double>& uR, vector<double>& uI, double phiR, double phiI)
{
	//int count = samples.size();
	//for (int i = 0; i < count; i++)
	//{
	//	sys->CalculateWavefunction(samples[i], uR, uI, phiR, phiI);
	//	for (int n = 0; n < MC_NTHERMSTEPS; n++)
	//	{
	//		DoMetropolisStep(samples[i], uR, uI, phiR, phiI);
	//	}
	//	wfValues2[i] = sys->GetWf() * sys->GetWf();
	//	exponentValues[i] = sys->GetExponent();
	//}
}

void UpdateSamples(vector<ICorrelatedSamplingData*>& samples, vector<double>& uR, vector<double>& uI, double phiR, double phiI)
{
	//int nrOfUpdatedSamples = 0;
	//int count = samples.size();
	//double weight = 1.0;
	//double weightDiff = 0.0;
	//
	//for (int i = 0; i < count; i++)
	//{
	//	sys->CalculateWavefunction(samples[i], uR, uI, phiR, phiI);
	//	//double newWf2 = sys->GetWf() * sys->GetWf();
	//	//weight = newWf2 / wfValues2[i];
	//	weight = exp(2.0 * sys->GetExponent() - 2.0 * correlatedSamplingData[i].exponent);
	//	weightDiff = abs(1.0 - weight);
	//
	//	if (weightDiff > 0.1)
	//	{
	//		nrOfUpdatedSamples++;
	//		for (int n = 0; n < MC_NTHERMSTEPS; n++)
	//		{
	//			DoMetropolisStep(samples[i], uR, uI, phiR, phiI);
	//		}
	//		wfValues2[i] = sys->GetWf() * sys->GetWf();
	//		exponentValues[i] = sys->GetExponent();
	//	}
	//}
	////Log("nrOfUpdatedSamples: " + to_string(nrOfUpdatedSamples) + "(" + to_string((double)nrOfUpdatedSamples / count * 100.0) + "%)");
	//double avgNrOfUpdatedSamples = MPIMethods::ReduceToAverage(&nrOfUpdatedSamples);
	//if (processRank == rootRank)
	//{
	//	Log("nrOfUpdatedSamples: " + to_string(avgNrOfUpdatedSamples) + "(" + to_string(avgNrOfUpdatedSamples / count * 100.0) + "%)");
	//}
}

/////////////////////////////////
/// update expectation values ///
/////////////////////////////////

void UpdateExpectationValues(vector<vector<double> >& R, vector<double>& uR, vector<double>& uI, double phiR, double phiI, bool intermediateStep = false)
{
	//vector<vector<double> > singlelocalOperators; // for all O_k
	vector<double> singlelocalEnergyR; // for all E^R
	//vector<double> singlelocalEnergyI; // for all E^I
	//vector<vector<vector<double> > > singlelocalOperatorsMatrix; // for all O_k O_j
	//vector<vector<double> > singlelocalOperatorlocalEnergyR; // for all O_k E^R
	//vector<vector<double> > singlelocalOperatorlocalEnergyI; // for all O_k E^I

	MC_NSTEPS *= (sys->GetStep() % WRITE_EVERY_NTH_STEP_TO_FILE == 0 ? MC_NSTEP_MULTIPLICATION_FACTOR_FOR_WRITE_DATA : 1);
	mc_nsteps = (double) MC_NSTEPS;
	ClearVector(localOperators);
	localEnergyR = 0;
	localEnergyI = 0;
	for (auto &row : localOperatorsMatrix)
	{
		ClearVector(row);
	}
	ClearVector(localOperatorlocalEnergyR);
	ClearVector(localOperatorlocalEnergyI);
	ClearVector(otherExpectationValues);

	sys->CalculateWavefunction(R, uR, uI, phiR, phiI);

	if (processRank == rootRank && !intermediateStep)
	{
		cout << "exponent=" << sys->GetExponent() << "\t\twf=" << sys->GetWf() << "\t\tphiR=" << phiR << endl;
	}
	//cout << "wf=" << Sys::wf << endl;
	//Sys::WriteLocalOperatorsToFile("start");
	for (int i = 0; i < MC_NINITIALIZATIONSTEPS; i++)
	{
		DoMetropolisStep(R, uR, uI, phiR, phiI);
	}
	//cout << "MC_NINITIALIZATIONSTEPS done" << endl;
	for (int i = 0; i < MC_NSTEPS; i++)
	{
		for (int nTherm = 0; nTherm < MC_NTHERMSTEPS; nTherm++)
		{
			DoMetropolisStep(R, uR, uI, phiR, phiI);
		}

		sys->CalculateExpectationValues(R, uR, uI, phiR, phiI);

		if (i < mc_nsteps_original)
		{
			correlatedSamplingData[i]->R = R;
			correlatedSamplingData[i]->wf = sys->GetWf();
			correlatedSamplingData[i]->wf2 = correlatedSamplingData[i]->wf * correlatedSamplingData[i]->wf;
			correlatedSamplingData[i]->exponent = sys->GetExponent();
			correlatedSamplingData[i]->exponent2 = correlatedSamplingData[i]->exponent * correlatedSamplingData[i]->exponent;
			correlatedSamplingData[i]->localOperators = sys->GetLocalOperators();
			sys->FillCorrelatedSamplingData(correlatedSamplingData[i]);
		}

		//INFO: consumes too much memory
		//singlelocalOperators.push_back(Sys::localOperators);
		singlelocalEnergyR.push_back(sys->GetLocalEnergyR());
		//singlelocalEnergyI.push_back(Sys::localEnergyI);
		//singlelocalOperatorsMatrix.push_back(Sys::localOperatorsMatrix);
		//singlelocalOperatorlocalEnergyR.push_back(Sys::localOperatorlocalEnergyR);
		//singlelocalOperatorlocalEnergyI.push_back(Sys::localOperatorlocalEnergyI);
		//WriteDataToFile(singlelocalEnergyR.back(), "singlelocalEnergyR", "er");

		//INFO: calculate contribution to average value
		localOperators += sys->GetLocalOperators() / mc_nsteps;
		localEnergyR += sys->GetLocalEnergyR() / mc_nsteps;
		localEnergyI += sys->GetLocalEnergyI() / mc_nsteps;
		localOperatorsMatrix += sys->GetLocalOperatorsMatrix() / mc_nsteps;
		localOperatorlocalEnergyR += sys->GetLocalOperatorlocalEnergyR() / mc_nsteps;
		localOperatorlocalEnergyI += sys->GetLocalOperatorlocalEnergyI() / mc_nsteps;
		otherExpectationValues += sys->GetOtherExpectationValues() / mc_nsteps;

		if ((10 * i) % MC_NSTEPS == 0)
		{
			//cout << (i / (double)MC_NSTEPS * 100.0) << "%" << endl;
			if (processRank == rootRank && !intermediateStep)
			{
				cout << "." << flush;
			}
		}
	}
	if (processRank == rootRank && !intermediateStep)
	{
		cout << endl;
	}

	//INFO: contribution to average values is calculated in each evaluation
	//localOperators = Mean(singlelocalOperators);
	//localEnergyR = Mean(singlelocalEnergyR);
	//localEnergyI = Mean(singlelocalEnergyI);
	//localOperatorsMatrix = Mean(singlelocalOperatorsMatrix);
	//localOperatorlocalEnergyR = Mean(singlelocalOperatorlocalEnergyR);
	//localOperatorlocalEnergyI = Mean(singlelocalOperatorlocalEnergyI);
	if (processRank == rootRank && !intermediateStep)
	{
		//WriteDataToFile(singlelocalEnergyR, "singlelocalEnergyR" + to_string(t), "singlelocalEnergyR");
		//Sys::WriteLocalOperatorsToFile(to_string(t));
		//Sys::WriteSplineSumsToFiles(to_string(t));
		//WriteParticlesToFile(R, to_string(t));
	}

	if (processRank == rootRank && !intermediateStep)
	{
		cout << "Acceptance: " << (nAcceptances / (nTrials / 100.0)) << "% (" << nAcceptances << "/" << nTrials << ")" << endl;
	}

	MC_NSTEPS /= (sys->GetStep() % WRITE_EVERY_NTH_STEP_TO_FILE == 0 ? MC_NSTEP_MULTIPLICATION_FACTOR_FOR_WRITE_DATA : 1);
}

void ParallelUpdateExpectationValues(vector<vector<double> >& R, vector<double>& uR, vector<double>& uI, double phiR, double phiI, bool intermediateStep = false)
{
	Timer t;
	double dblDuration;
	t.start();

	UpdateExpectationValues(R, uR, uI, phiR, phiI, intermediateStep);

	t.stop();
	dblDuration = (double) t.duration();
	vector<double> timings = MPIMethods::ReduceToMinMaxMean(dblDuration);
	if (processRank == rootRank && !intermediateStep)
	{
		Log("duration: min = " + to_string(timings[0]) + " ms");
		Log("          max = " + to_string(timings[1]) + " ms");
		Log("          <t> = " + to_string(timings[2]) + " ms");
	}

	//Log("localOperators" + to_string(processRank));
	//WriteDataToFile(localOperators, "localOperators" + to_string(processRank), "localOperators");

	MPIMethods::ReduceToAverage(localOperators);
	MPIMethods::ReduceToAverage(&localEnergyR);
	MPIMethods::ReduceToAverage(&localEnergyI);
	MPIMethods::ReduceToAverage(localOperatorsMatrix);
	MPIMethods::ReduceToAverage(localOperatorlocalEnergyR);
	MPIMethods::ReduceToAverage(localOperatorlocalEnergyI);
	MPIMethods::ReduceToAverage(otherExpectationValues);

	//BulkOnlySplines* s = dynamic_cast<BulkOnlySplines*>(sys);
	//vector<vector<double> > localKineticEnergiesD1 = Mean(s->allLocalKineticEnergiesD1);
	//vector<double> localKineticEnergiesD2 = Mean(s->allLocalKineticEnergiesD2);
	//if (processRank == rootRank)
	//{
	//	WriteDataToFile(s->allER1, "BB_allER1", to_string(Mean(s->allER1)));
	//	WriteDataToFile(s->allER2, "BB_allER2", to_string(Mean(s->allER2)));
	//	WriteDataToFile(s->allER1new, "BB_allER1new", to_string(Mean(s->allER1new)));
	//	WriteDataToFile(s->allER2new, "BB_allER2new", to_string(Mean(s->allER2new)));
	//	WriteDataToFile(s->splineSumsD, "BB_splineSumsD", "test");
	//	WriteDataToFile(s->splineSumsD2, "BB_splineSumsD2", "test");
	//	WriteDataToFile(s->allLocalKineticEnergiesD1, "BB_Original_localKineticEnergiesD1", "test");
	//	WriteDataToFile(s->allLocalKineticEnergiesD2, "BB_Original_localKineticEnergiesD2", "test");
	//	WriteDataToFile(localKineticEnergiesD1, "BB_Mean_localKineticEnergiesD1", "test", 1);
	//	WriteDataToFile(localKineticEnergiesD2, "BB_Mean_localKineticEnergiesD2", to_string(Sum(localKineticEnergiesD2)), 1);
	//	auto t = OuterSum(s->allLocalKineticEnergiesD2);
	//	WriteDataToFile(t, "BB_OuterSum_allLocalKineticEnergiesD2", "outersum", 1);
	//	auto g = OuterSum(s->splineSumsOnlyD2);
	//	WriteDataToFile(g, "BB_splineSumsOnlyD2", to_string(Sum(g)), 1);
	//}
	//auto os2 = OuterSum(s->allLocalKineticEnergiesD2);
	//MPIMethods::ReduceToAverage(os2);
	////MPIMethods::ReduceToAverage(localKineticEnergiesD1);
	////MPIMethods::ReduceToAverage(localKineticEnergiesD2);
	//if (processRank == rootRank)
	//{
	//	WriteDataToFile(localKineticEnergiesD1, "BB_localKineticEnergiesD1", "test", 1);
	//	WriteDataToFile(localKineticEnergiesD2, "BB_localKineticEnergiesD2", to_string(Sum(localKineticEnergiesD2)), 1);
	//	WriteDataToFile(os2, "BB_OuterSum_allLocalKineticEnergiesD2_Reduced", to_string(Sum(os2)), 1);
	//}
}

void UpdateExpectationValuesForGivenSamples(vector<ICorrelatedSamplingData*>& samples, vector<double>& uR, vector<double>& uI, double phiR, double phiI, bool intermediateStep = false)
{
	int count = samples.size();
	double dblCount = (double) count;
	//double weight = 1.0;
	//double weightSum = 0.0;

	//vector<double> newWf2s;
	//vector<double> weights;
	//vector<vector<double> > los;

	ClearVector(localOperators);
	localEnergyR = 0;
	localEnergyI = 0;
	for (auto &row : localOperatorsMatrix)
	{
		ClearVector(row);
	}
	ClearVector(localOperatorlocalEnergyR);
	ClearVector(localOperatorlocalEnergyI);
	ClearVector(otherExpectationValues);

	for (int i = 0; i < count; i++)
	{
		sys->CalculateWavefunction(samples[i], uR, uI, phiR, phiI);
		sys->CalculateExpectationValues(samples[i], uR, uI, phiR, phiI);

		//weight = sys->GetWf() * sys->GetWf() / wfValues2[i];
		//weight = exp(2.0 * sys->GetExponent() - 2.0 * samples[i]->exponent);
		//weight = sys->GetExponent() / samples[i]->exponent;
		//weight = 1.0;
		//weights.push_back(weight);
		//weightSum += weight;
		//INFO: calculate contribution to average value
		localOperators += sys->GetLocalOperators() / dblCount;		// * weight;
		localEnergyR += sys->GetLocalEnergyR() / dblCount;		// * weight;
		localEnergyI += sys->GetLocalEnergyI() / dblCount;		// * weight;
		localOperatorsMatrix += sys->GetLocalOperatorsMatrix() / dblCount;		// * weight;
		localOperatorlocalEnergyR += sys->GetLocalOperatorlocalEnergyR() / dblCount;		// * weight;
		localOperatorlocalEnergyI += sys->GetLocalOperatorlocalEnergyI() / dblCount;		// * weight;
		otherExpectationValues += sys->GetOtherExpectationValues() / dblCount;		// * weight;

		//weights.push_back(weight);
		//newWf2s.push_back(sys->GetWf() * sys->GetWf());
		//los.push_back(sys->GetLocalOperators());

		if ((10 * i) % count == 0)
		{
			//cout << (i / (double)MC_NSTEPS * 100.0) << "%" << endl;
			if (processRank == rootRank && !intermediateStep)
			{
				cout << "." << flush;
				//cout << weight << endl;
			}
		}
	}

	//if (processRank == rootRank && sys->GetStep() > 1 && !intermediateStep)
	//{
	//	WriteDataToFile(weights, "__Test_weights", "phiR=" + to_string(phiR));
	//	WriteDataToFile(los, "__Test_los", "los");
	//	WriteDataToFile(newWf2s, "__Test_newWf2s", "newWf2s");
	//}

	//double inverseWeightMean = dblCount / weightSum;
	//if (processRank == rootRank && !intermediateStep)
	//{
	//	cout << "mean: " << (1.0/inverseWeightMean) << endl;
	//}
	//localOperators *= inverseWeightMean;
	//localEnergyR *= inverseWeightMean;
	//localEnergyI *= inverseWeightMean;
	//localOperatorsMatrix *= inverseWeightMean;
	//localOperatorlocalEnergyR *= inverseWeightMean;
	//localOperatorlocalEnergyI *= inverseWeightMean;
	//otherExpectationValues *= inverseWeightMean;

	if (processRank == rootRank && !intermediateStep)
	{
		cout << endl;
	}
}

void ParallelUpdateExpectationValuesForGivenSamples(vector<ICorrelatedSamplingData*>& samples, vector<double>& uR, vector<double>& uI, double phiR, double phiI, bool intermediateStep = false)
{
	//Timer t;
	//double dblDuration;
	//t.start();

	UpdateExpectationValuesForGivenSamples(samples, uR, uI, phiR, phiI, intermediateStep);

	//t.stop();
	//dblDuration = (double) t.duration();
	//vector<double> timings = MPIMethods::ReduceToMinMaxMean(dblDuration);
	//if (processRank == rootRank)
	//{
	//	Log("duration: min = " + to_string(timings[0]) + " ms");
	//	Log("          max = " + to_string(timings[1]) + " ms");
	//	Log("          <t> = " + to_string(timings[2]) + " ms");
	//}

	MPIMethods::ReduceToAverage(localOperators);
	MPIMethods::ReduceToAverage(&localEnergyR);
	MPIMethods::ReduceToAverage(&localEnergyI);
	MPIMethods::ReduceToAverage(localOperatorsMatrix);
	MPIMethods::ReduceToAverage(localOperatorlocalEnergyR);
	MPIMethods::ReduceToAverage(localOperatorlocalEnergyI);
	MPIMethods::ReduceToAverage(otherExpectationValues);
}

void CalculateAdditionalSystemProperties(vector<vector<double> >& R, vector<double>& uR, vector<double>& uI, double phiR, double phiI)
{
	int percent = 0;
	ClearVector(additionalSystemProperties);
	ClearVector(AllAdditionalSystemProperties);
	sys->CalculateWavefunction(R, uR, uI, phiR, phiI);
	for (int i = 0; i < MC_NADDITIONALINITIALIZATIONSTEPS; i++)
	{
		DoMetropolisStep(R, uR, uI, phiR, phiI);
	}
	for (int i = 0; i < MC_NADDITIONALSTEPS; i++)
	{
		for (int nTherm = 0; nTherm < MC_NADDITIONALTHERMSTEPS; nTherm++)
		{
			DoMetropolisStep(R, uR, uI, phiR, phiI);
		}

		sys->CalculateAdditionalSystemProperties(R, uR, uI, phiR, phiI);

		//INFO: calculate contribution to average value
		additionalSystemProperties += sys->GetAdditionalSystemProperties() / mc_nadditionalsteps;
		AllAdditionalSystemProperties.push_back(sys->GetAdditionalSystemProperties());

		if ((100 * i) % MC_NADDITIONALSTEPS == 0)
		{
			//cout << (i / (double)MC_NSTEPS * 100.0) << "%" << endl;
			if (processRank == rootRank)
			{
				//cout << "." << flush;
				percent++;
				cout << percent << "%" << endl;
			}
		}
	}
	if (processRank == rootRank)
	{
		cout << endl;
	}
	if (processRank == rootRank)
	{
		cout << "Acceptance: " << (nAcceptances / (nTrials / 100.0)) << "% (" << nAcceptances << "/" << nTrials << ")" << endl;
	}
}

void ParallelCalculateAdditionalSystemProperties(vector<vector<double> >& R, vector<double>& uR, vector<double>& uI, double phiR, double phiI)
{
	CalculateAdditionalSystemProperties(R, uR, uI, phiR, phiI);

	MPIMethods::ReduceToAverage(additionalSystemProperties);
}

/////////////////////////////////
/// build system of equations ///
/////////////////////////////////

void BuildSystemOfEquationsForParametersNoPhi(vector<vector<double> >& matrix, vector<double>& energiesReal, vector<double>& energiesImag)
{
	energiesReal.resize(N_PARAM);
	energiesImag.resize(N_PARAM);
	matrix.resize(N_PARAM);
	for (auto &i : matrix)
	{
		i.resize(N_PARAM);
	}

	for (int i = 0; i < N_PARAM; i++)
	{
		if (IMAGINARY_TIME == 0)
		{
			energiesReal[i] = localOperatorlocalEnergyI[i];
			energiesImag[i] = -localOperatorlocalEnergyR[i];
		}
		else
		{
			energiesReal[i] = -localOperatorlocalEnergyR[i];
			energiesImag[i] = -localOperatorlocalEnergyI[i];
		}
		for (int j = 0; j <= i; j++)
		{
			matrix[i][j] = localOperatorsMatrix[i][j];
			matrix[j][i] = matrix[i][j];
		}
	}
}

void BuildSystemOfEquationsForParametersIncludePhiWithTimeRotation(vector<vector<double> >& matrix, vector<double>& energiesReal, vector<double>& energiesImag)
{
	double rotation = 1.499 * M_PI; // 3/2 Pi -> real time; Pi -> imaginary time
	double cosRotation = cos(rotation);
	double sinRotation = sin(rotation);

	energiesReal.resize(N_PARAM);
	energiesImag.resize(N_PARAM);
	matrix.resize(N_PARAM);
	for (auto &i : matrix)
	{
		i.resize(N_PARAM);
	}

	for (int i = 0; i < N_PARAM; i++)
	{
		energiesReal[i] = cosRotation * (localOperatorlocalEnergyR[i] - localEnergyR * localOperators[i]) - sinRotation * (localOperatorlocalEnergyI[i]);
		energiesImag[i] = sinRotation * (localOperatorlocalEnergyR[i] - localEnergyR * localOperators[i]) + cosRotation * (localOperatorlocalEnergyI[i]);
		for (int j = 0; j <= i; j++)
		{
			matrix[i][j] = localOperatorsMatrix[i][j] - localOperators[i] * localOperators[j];
			matrix[j][i] = matrix[i][j];
		}
	}
}

void BuildSystemOfEquationsForParametersIncludePhi(vector<vector<double> >& matrix, vector<double>& energiesReal, vector<double>& energiesImag)
{
	energiesReal.resize(N_PARAM);
	energiesImag.resize(N_PARAM);
	matrix.resize(N_PARAM);
	for (auto &i : matrix)
	{
		i.resize(N_PARAM);
	}

	for (int i = 0; i < N_PARAM; i++)
	{
		if (IMAGINARY_TIME == 0)
		{
			energiesReal[i] = localOperatorlocalEnergyI[i];
			energiesImag[i] = -localOperatorlocalEnergyR[i] + localEnergyR * localOperators[i];
		}
		else
		{
			energiesReal[i] = -localOperatorlocalEnergyR[i] + localEnergyR * localOperators[i];
			energiesImag[i] = -localOperatorlocalEnergyI[i];
		}
		for (int j = 0; j <= i; j++)
		{
			matrix[i][j] = localOperatorsMatrix[i][j] - localOperators[i] * localOperators[j];
			matrix[j][i] = matrix[i][j];
		}
	}
}

void BuildSystemOfEquationsForParameters(vector<vector<double> >& matrix, vector<double>& energiesReal, vector<double>& energiesImag)
{
	if (IMAGINARY_TIME == -1)
	{
		BuildSystemOfEquationsForParametersIncludePhiWithTimeRotation(matrix, energiesReal, energiesImag);
	}
	else if (sys->USE_NORMALIZATION_AND_PHASE)
	{
		BuildSystemOfEquationsForParametersIncludePhi(matrix, energiesReal, energiesImag);
	}
	else
	{
		//BuildSystemOfEquationsForParametersNoPhi(matrix, energiesReal, energiesImag);
		BuildSystemOfEquationsForParametersIncludePhi(matrix, energiesReal, energiesImag);
	}
}

/////////////////////////////////
/// solve system of equations ///
/////////////////////////////////

void PerformCholeskyDecomposition(vector<vector<double> >& matrix)
{
	double sum = 0;
	bool logged = false;
	for (int i = 0; i < N_PARAM; i++)
	{
		for (int j = 0; j <= i; j++)
		{
			sum = matrix[i][j];
			for (int k = 0; k < j; k++)
			{
				sum -= matrix[i][k] * matrix[j][k];
			}
			if (i > j)
			{
				matrix[i][j] = sum / matrix[j][j];
			}
			else if (sum > 0)
			{
				matrix[i][i] = sqrt(sum);
			}
			else
			{
				if (!logged)
				{
					Log("!!!!! NOT POSITIVE SEMI DEFINITE !!!! @ i=" + to_string(i) + ", j=" + to_string(j), ERROR);
					doNotAcceptStep = true;
					logged = true;
				}
			}
		}
	}
}

void SolveCholeskyDecomposedEquationSystem(vector<vector<double> >& matrix, vector<double>& rhs, vector<double>& solution)
{
	int size = rhs.size();
	double sum = 0;
	vector<double> tmp;
	tmp.resize(size);

	solution.resize(size);

	for (int i = 0; i < N_PARAM; i++)
	{
		sum = 0;
		for (int j = 0; j < i; j++)
		{
			sum += matrix[i][j] * tmp[j];
		}
		tmp[i] = 1.0 / matrix[i][i] * (rhs[i] - sum);
	}

	for (int i = N_PARAM - 1; i >= 0; i--)
	{
		sum = 0;
		for (int j = N_PARAM - 1; j > i; j--)
		{
			sum += matrix[j][i] * solution[j];
		}
		solution[i] = 1.0 / matrix[i][i] * (tmp[i] - sum);
	}
}

void SolveEquationSystem(vector<vector<double> >& matrix, vector<double>& rhs, vector<double>& solution)
{
	PerformCholeskyDecomposition(matrix);
	SolveCholeskyDecomposedEquationSystem(matrix, rhs, solution);
}

void SolveEquationSystemArmadillo(vector<vector<double> >& matrix, vector<double>& rhs, vector<double>& solution)
{
	cout << "#############################################" << endl;
	cout << "#############################################" << endl;
	cout << "SolveEquationSystemArmadillo NOT IMPLEMENTED!" << endl;
	cout << "#############################################" << endl;
	cout << "#############################################" << endl;
	//arma::mat m(N_PARAM, N_PARAM);
	//arma::vec b(N_PARAM);
	//arma::vec x(N_PARAM);
	//
	//for (int i = 0; i < N_PARAM; i++)
	//{
	//	for (int j = 0; j < N_PARAM; j++)
	//	{
	//		m(i, j) = matrix[i][j];
	//	}
	//	b[i] = rhs[i];
	//}
	//
	//x = arma::solve(m, b);
	//
	//for (int i = 0; i < N_PARAM; i++)
	//{
	//	solution[i] = x[i];
	//}
}

void CalculatePhiDot(vector<double>& uDotR, vector<double>& uDotI, double *phiDotR, double *phiDotI)
{
	for (int i = 0; i < N_PARAM; i++)
	{
		*phiDotR -= localOperators[i] * uDotR[i];
		*phiDotI -= localOperators[i] * uDotI[i];
	}

	if (IMAGINARY_TIME == -1)
	{
		double rotation = 1.499 * M_PI; // 3/2 Pi -> real time; Pi -> imaginary time
		double cosRotation = cos(rotation);
		double sinRotation = sin(rotation);
		*phiDotI -= cosRotation * localEnergyR;
		*phiDotR -= sinRotation * localEnergyR;
	}
	else if (IMAGINARY_TIME == 0)
	{
		*phiDotI -= localEnergyR;
	}
	else //INFO: IMAGINARY_TIME == 1
	{
		*phiDotR -= localEnergyR;
	}
}

void PreconditionEquationSystemByScaling(vector<vector<double> >& matrix, vector<double>& energiesReal, vector<double>& energiesImag, vector<double>& preconditionScalings)
{
	//INFO: for square matrices only!
	preconditionScalings.resize(matrix.size());
	for (unsigned int i = 0; i < matrix.size(); i++)
	{
		preconditionScalings[i] = sqrt(matrix[i][i]);
	}
	for (unsigned int i = 0; i < matrix.size(); i++)
	{
		for (unsigned int j = 0; j < matrix.size(); j++)
		{
			matrix[i][j] /= preconditionScalings[i] * preconditionScalings[j];
		}
		energiesReal[i] /= preconditionScalings[i];
		energiesImag[i] /= preconditionScalings[i];
	}
}

void RegularizeEquationSystem(vector<vector<double> >& matrix, double eps)
{
	//INFO: for square matrices only!
	for (unsigned int i = 0; i < matrix.size(); i++)
	{
		matrix[i][i] += eps;
	}
}

void SolveForParametersDot(vector<double>& uDotR, vector<double>& uDotI, double *phiDotR, double *phiDotI)
{
	vector<double> energiesReal; //rhs of the equation that corresponds to uDotR;
	vector<double> energiesImag; //rhs of the equation that corresponds to uDotI;
	vector<vector<double> > matrix;

	BuildSystemOfEquationsForParameters(matrix, energiesReal, energiesImag);

	//WriteDataToFile(matrix, "Omatrix", "Omatrix");

	vector<double> preconditionScalings;
	if (USE_PRECONDITIONING == 1)
	{
		PreconditionEquationSystemByScaling(matrix, energiesReal, energiesImag, preconditionScalings);
		RegularizeEquationSystem(matrix, 0.001);
	}

	//WriteDataToFile(matrix, "PCmatrix", "PCmatrix");
	//WriteDataToFile(energiesReal, "BBenergiesReal", "BBenergiesReal");
	//WriteDataToFile(energiesImag, "BBenergiesImag", "BBenergiesImag");

	PerformCholeskyDecomposition(matrix);
	SolveCholeskyDecomposedEquationSystem(matrix, energiesReal, uDotR);
	SolveCholeskyDecomposedEquationSystem(matrix, energiesImag, uDotI);
	CalculatePhiDot(uDotR, uDotI, phiDotR, phiDotI);

	if (USE_PRECONDITIONING)
	{
		for (unsigned int i = 0; i < uDotR.size(); i++)
		{
			uDotR[i] /= preconditionScalings[i];
			uDotI[i] /= preconditionScalings[i];
		}
	}

	//WriteDataToFile(matrix, "BBcholesky", "BBcholesky");
	//WriteDataToFile(uDotR, "BBuDotR", "BBuDotR");
}

///////////////////////
/// ODE integrators ///
///////////////////////

void CalculateNextParametersEuler(double dt, vector<double>& uR, vector<double>& uI, double *phiR, double *phiI)
{
	if (processRank == rootRank)
	{
		vector<double> uDotR;
		vector<double> uDotI;
		double phiDotR = 0;
		double phiDotI = 0;

		SolveForParametersDot(uDotR, uDotI, &phiDotR, &phiDotI);

		uR = uR + uDotR * dt;
		uI = uI + uDotI * dt;
		*phiR = *phiR + phiDotR * dt;
		*phiI = *phiI + phiDotI * dt;

		//uR[N_PARAM - 2] += uDotR[N_PARAM - 2] * 9.0 * dt;
		//uR[N_PARAM - 1] += uDotR[N_PARAM - 1] * 4.0 * dt;
	}
}

void CalculateNextParametersPC(double dt, vector<vector<double> >& R, vector<double>& uR, vector<double>& uI, double *phiR, double *phiI)
{
	double dt_2 = dt / 2.0;

	vector<double> uDotR;
	vector<double> uDotI;
	double phiDotR = 0;
	double phiDotI = 0;

	int PCsteps = 1;
	vector<double> tmpUR(N_PARAM);
	vector<double> tmpUI(N_PARAM);
	double tmpPhiR = 0;
	double tmpPhiI = 0;
	vector<double> nextUDotR;
	vector<double> nextUDotI;
	double nextPhiDotR = 0;
	double nextPhiDotI = 0;

	ClearVector(tmpUR);
	ClearVector(tmpUI);

	if (processRank == rootRank)
	{
		SolveForParametersDot(uDotR, uDotI, &phiDotR, &phiDotI);

		tmpUR = uR + uDotR * dt;
		tmpUI = uI + uDotI * dt;
		tmpPhiR = *phiR + phiDotR * dt;
		tmpPhiI = *phiI + phiDotI * dt;
	}
	for (int s = 0; s < PCsteps; s++)
	{
		//BroadcastNewParameters(uR, uI, phiR, phiI);
		//BroadcastValues(tmpUR, N_PARAM);
		//BroadcastValues(tmpUI, N_PARAM);
		//BroadcastValue(&tmpPhiR);
		//BroadcastValue(&tmpPhiI);
		BroadcastNewParameters(tmpUR, tmpUI, &tmpPhiR, &tmpPhiI);
		ParallelUpdateExpectationValues(R, tmpUR, tmpUI, tmpPhiR, tmpPhiI, true);
		if (processRank == rootRank)
		{
			SolveForParametersDot(nextUDotR, nextUDotI, &nextPhiDotR, &nextPhiDotI);

			tmpUR = uR + (uDotR + nextUDotR) * dt_2;
			tmpUI = uI + (uDotI + nextUDotI) * dt_2;
			tmpPhiR = *phiR + (phiDotR + nextPhiDotR) * dt_2;
			tmpPhiI = *phiI + (phiDotI + nextPhiDotI) * dt_2;
		}
	}

	uR = tmpUR;
	uI = tmpUI;
	*phiR = tmpPhiR;
	*phiI = tmpPhiI;
}

void CalculateNextParametersPCReuseSamples(double dt, vector<double>& uR, vector<double>& uI, double *phiR, double *phiI)
{
	double dt_2 = dt / 2.0;

	vector<double> uDotR;
	vector<double> uDotI;
	double phiDotR = 0;
	double phiDotI = 0;

	int PCsteps = 6;
	vector<double> tmpUR(N_PARAM);
	vector<double> tmpUI(N_PARAM);
	double tmpPhiR = 0;
	double tmpPhiI = 0;
	vector<double> nextUDotR;
	vector<double> nextUDotI;
	double nextPhiDotR = 0;
	double nextPhiDotI = 0;

	ClearVector(tmpUR);
	ClearVector(tmpUI);

	if (processRank == rootRank)
	{
		SolveForParametersDot(uDotR, uDotI, &phiDotR, &phiDotI);

		tmpUR = uR + uDotR * dt;
		tmpUI = uI + uDotI * dt;
		tmpPhiR = *phiR + phiDotR * dt;
		tmpPhiI = *phiI + phiDotI * dt;
	}
	for (int s = 0; s < PCsteps; s++)
	{
		//BroadcastNewParameters(uR, uI, phiR, phiI);
		//BroadcastValues(tmpUR, N_PARAM);
		//BroadcastValues(tmpUI, N_PARAM);
		//BroadcastValue(&tmpPhiR);
		//BroadcastValue(&tmpPhiI);
		BroadcastNewParameters(tmpUR, tmpUI, &tmpPhiR, &tmpPhiI);
		ParallelUpdateExpectationValuesForGivenSamples(correlatedSamplingData, tmpUR, tmpUI, tmpPhiR, tmpPhiI, true);
		if (processRank == rootRank)
		{
			SolveForParametersDot(nextUDotR, nextUDotI, &nextPhiDotR, &nextPhiDotI);

			tmpUR = uR + (uDotR + nextUDotR) * dt_2;
			tmpUI = uI + (uDotI + nextUDotI) * dt_2;
			tmpPhiR = *phiR + (phiDotR + nextPhiDotR) * dt_2;
			tmpPhiI = *phiI + (phiDotI + nextPhiDotI) * dt_2;
		}
	}

	uR = tmpUR;
	uI = tmpUI;
	*phiR = tmpPhiR;
	*phiI = tmpPhiI;
}

void CalculateNextParametersRK4(double dt, vector<vector<double> >& R, vector<double>& uR, vector<double>& uI, double *phiR, double *phiI)
{
	double dt_2 = dt / 2.0;

	vector<double> tmpUR(N_PARAM);
	vector<double> tmpUI(N_PARAM);
	double tmpPhiR = 0;
	double tmpPhiI = 0;

	vector<vector<double> > uDotR;
	vector<vector<double> > uDotI;
	vector<double> phiDotR;
	vector<double> phiDotI;

	ClearVector(tmpUR);
	ClearVector(tmpUI);
	uDotR.resize(4);
	uDotI.resize(4);
	phiDotR.resize(4);
	phiDotI.resize(4);

	if (processRank == rootRank)
	{
		SolveForParametersDot(uDotR[0], uDotI[0], &(phiDotR[0]), &(phiDotI[0]));

		tmpUR = uR + uDotR[0] * dt_2;
		tmpUI = uI + uDotI[0] * dt_2;
		tmpPhiR = *phiR + phiDotR[0] * dt_2;
		tmpPhiI = *phiI + phiDotI[0] * dt_2;
	}
	BroadcastNewParameters(tmpUR, tmpUI, &tmpPhiR, &tmpPhiI);
	ParallelUpdateExpectationValues(R, tmpUR, tmpUI, tmpPhiR, tmpPhiI, true);

	if (processRank == rootRank)
	{
		SolveForParametersDot(uDotR[1], uDotI[1], &(phiDotR[1]), &(phiDotI[1]));

		tmpUR = uR + uDotR[1] * dt_2;
		tmpUI = uI + uDotI[1] * dt_2;
		tmpPhiR = *phiR + phiDotR[1] * dt_2;
		tmpPhiI = *phiI + phiDotI[1] * dt_2;
	}
	BroadcastNewParameters(tmpUR, tmpUI, &tmpPhiR, &tmpPhiI);
	ParallelUpdateExpectationValues(R, tmpUR, tmpUI, tmpPhiR, tmpPhiI, true);

	if (processRank == rootRank)
	{
		SolveForParametersDot(uDotR[2], uDotI[2], &(phiDotR[2]), &(phiDotI[2]));

		tmpUR = uR + uDotR[2] * dt;
		tmpUI = uI + uDotI[2] * dt;
		tmpPhiR = *phiR + phiDotR[2] * dt;
		tmpPhiI = *phiI + phiDotI[2] * dt;
	}
	BroadcastNewParameters(tmpUR, tmpUI, &tmpPhiR, &tmpPhiI);
	ParallelUpdateExpectationValues(R, tmpUR, tmpUI, tmpPhiR, tmpPhiI, true);

	if (processRank == rootRank)
	{
		SolveForParametersDot(uDotR[3], uDotI[3], &(phiDotR[3]), &(phiDotI[3]));

		uR = uR + ((uDotR[0] + uDotR[1] * 2.0 + uDotR[2] * 2.0 + uDotR[3]) / 6.0) * dt;
		uI = uI + ((uDotI[0] + uDotI[1] * 2.0 + uDotI[2] * 2.0 + uDotI[3]) / 6.0) * dt;
		*phiR = *phiR + ((phiDotR[0] + phiDotR[1] * 2.0 + phiDotR[2] * 2.0 + phiDotR[3]) / 6.0) * dt;
		*phiI = *phiI + ((phiDotI[0] + phiDotI[1] * 2.0 + phiDotI[2] * 2.0 + phiDotI[3]) / 6.0) * dt;
	}
}

void CalculateNextParametersRK4ReuseSamples(double dt, vector<double>& uR, vector<double>& uI, double *phiR, double *phiI)
{
	double dt_2 = dt / 2.0;

	vector<double> tmpUR(N_PARAM);
	vector<double> tmpUI(N_PARAM);
	double tmpPhiR = 0;
	double tmpPhiI = 0;

	vector<vector<double> > uDotR;
	vector<vector<double> > uDotI;
	vector<double> phiDotR;
	vector<double> phiDotI;

	ClearVector(tmpUR);
	ClearVector(tmpUI);
	uDotR.resize(4);
	uDotI.resize(4);
	phiDotR.resize(4);
	phiDotI.resize(4);

	if (processRank == rootRank)
	{
		SolveForParametersDot(uDotR[0], uDotI[0], &(phiDotR[0]), &(phiDotI[0]));

		tmpUR = uR + uDotR[0] * dt_2;
		tmpUI = uI + uDotI[0] * dt_2;
		tmpPhiR = *phiR + phiDotR[0] * dt_2;
		tmpPhiI = *phiI + phiDotI[0] * dt_2;
	}
	BroadcastNewParameters(tmpUR, tmpUI, &tmpPhiR, &tmpPhiI);
	ParallelUpdateExpectationValuesForGivenSamples(correlatedSamplingData, tmpUR, tmpUI, tmpPhiR, tmpPhiI, true);

	if (processRank == rootRank)
	{
		SolveForParametersDot(uDotR[1], uDotI[1], &(phiDotR[1]), &(phiDotI[1]));

		tmpUR = uR + uDotR[1] * dt_2;
		tmpUI = uI + uDotI[1] * dt_2;
		tmpPhiR = *phiR + phiDotR[1] * dt_2;
		tmpPhiI = *phiI + phiDotI[1] * dt_2;
	}
	BroadcastNewParameters(tmpUR, tmpUI, &tmpPhiR, &tmpPhiI);
	ParallelUpdateExpectationValuesForGivenSamples(correlatedSamplingData, tmpUR, tmpUI, tmpPhiR, tmpPhiI, true);

	if (processRank == rootRank)
	{
		SolveForParametersDot(uDotR[2], uDotI[2], &(phiDotR[2]), &(phiDotI[2]));

		tmpUR = uR + uDotR[2] * dt;
		tmpUI = uI + uDotI[2] * dt;
		tmpPhiR = *phiR + phiDotR[2] * dt;
		tmpPhiI = *phiI + phiDotI[2] * dt;
	}
	BroadcastNewParameters(tmpUR, tmpUI, &tmpPhiR, &tmpPhiI);
	ParallelUpdateExpectationValuesForGivenSamples(correlatedSamplingData, tmpUR, tmpUI, tmpPhiR, tmpPhiI, true);

	if (processRank == rootRank)
	{
		SolveForParametersDot(uDotR[3], uDotI[3], &(phiDotR[3]), &(phiDotI[3]));

		uR = uR + ((uDotR[0] + uDotR[1] * 2.0 + uDotR[2] * 2.0 + uDotR[3]) / 6.0) * dt;
		uI = uI + ((uDotI[0] + uDotI[1] * 2.0 + uDotI[2] * 2.0 + uDotI[3]) / 6.0) * dt;
		*phiR = *phiR + ((phiDotR[0] + 2.0 * phiDotR[1] + 2.0 * phiDotR[2] + phiDotR[3]) / 6.0) * dt;
		*phiI = *phiI + ((phiDotI[0] + 2.0 * phiDotI[1] + 2.0 * phiDotI[2] + phiDotI[3]) / 6.0) * dt;
	}
}

void CalculateNextParametersImplicitEuler(double dt, vector<double>& uR, vector<double>& uI, double *phiR, double *phiI)
{
	vector<double> uR_Start(uR);
	vector<double> uI_Start(uI);
	double phiR_Start = *phiR;
	double phiI_Start = *phiI;

	vector<double> uDotR(N_PARAM);
	vector<double> uDotI(N_PARAM);
	double phiDotR = 0;
	double phiDotI = 0;

	vector<double> uDotR_n(N_PARAM);
	vector<double> uDotI_n(N_PARAM);
	double phiDotR_n = 0;
	double phiDotI_n = 0;

	vector<double> uDotR_n_delta(N_PARAM);
	vector<double> uDotI_n_delta(N_PARAM);
	double phiDotR_n_delta = 0;
	double phiDotI_n_delta = 0;

	vector<vector<double> > JR(N_PARAM);
	vector<vector<double> > JI(N_PARAM);
	for (int i = 0; i < N_PARAM; i++)
	{
		JR[i].resize(N_PARAM);
		JI[i].resize(N_PARAM);
	}

	if (processRank == rootRank)
	{
		SolveForParametersDot(uDotR, uDotI, &phiDotR, &phiDotI);

		uR_Start = uR_Start + uDotR * dt;
		uI_Start = uI_Start + uDotI * dt;
		phiR_Start = phiR_Start + phiDotR * dt;
		phiI_Start = phiI_Start + phiDotI * dt;
	}

	BroadcastNewParameters(uR_Start, uI_Start, &phiR_Start, &phiI_Start);
	BroadcastNewParameters(uDotR, uDotI, &phiDotR, &phiDotI);

	double eps = 1e-2;
	double diff = 1;
	double delta = 1e-4;
	int iterations = 0;

	vector<double> uR_n(uR_Start);
	vector<double> uI_n(uI_Start);
	double phiR_n = phiR_Start;
	double phiI_n = phiI_Start;
	//vector<double> uR_n(uR);
	//vector<double> uI_n(uI);
	//double phiR_n = *phiR;
	//double phiI_n = *phiI;

	//vector<double> directionR(N_PARAM);
	//vector<double> directionI(N_PARAM);
	//for (int i = 0; i < N_PARAM; i++)
	//{
	//	directionR[i] = uDotR[i] > 0 ? 1 : -1;
	//	directionI[i] = uDotI[i] > 0 ? 1 : -1;
	//}

	while (diff > eps && iterations < 20)
	{
		iterations++;
		ParallelUpdateExpectationValuesForGivenSamples(correlatedSamplingData, uR_n, uI_n, phiR_n, phiI_n, true);
		if (processRank == rootRank)
		{
			SolveForParametersDot(uDotR_n, uDotI_n, &phiDotR_n, &phiDotI_n);
		}
		BroadcastNewParameters(uDotR_n, uDotI_n, &phiDotR_n, &phiDotI_n);

		vector<double> gR_n = uR_n - uR - uDotR_n * dt; //TODO: gR_n and uDotR_n is only used by rootRank?!?
		vector<double> gI_n = uI_n - uI - uDotI_n * dt;

		for (int i = 0; i < N_PARAM; i++)
		{
			double dUR = delta * uR[i]; // * directionR[i];
			uR_n[i] += dUR;
			ParallelUpdateExpectationValuesForGivenSamples(correlatedSamplingData, uR_n, uI_n, phiR_n, phiI_n, true);
			if (processRank == rootRank)
			{
				SolveForParametersDot(uDotR_n_delta, uDotI_n_delta, &phiDotR_n_delta, &phiDotI_n_delta);
				JR[i] = (uDotR_n_delta - uDotR_n) * dt / dUR * (-1.0);
				JR[i][i] = JR[i][i] + 1.0;
			}
			uR_n[i] -= dUR;

			double dUI = delta * uI[i]; // * directionI[i];
			uI_n[i] += dUI;
			ParallelUpdateExpectationValuesForGivenSamples(correlatedSamplingData, uR_n, uI_n, phiR_n, phiI_n, true);
			if (processRank == rootRank)
			{
				SolveForParametersDot(uDotR_n_delta, uDotI_n_delta, &phiDotR_n_delta, &phiDotI_n_delta);
				JI[i] = (uDotI_n_delta - uDotI_n) * dt / dUI * (-1.0);
				JI[i][i] = JI[i][i] + 1.0;
			}
			uI_n[i] -= dUI;
		}
		if (processRank == rootRank)
		{
			vector<double> delta_uR(N_PARAM);
			vector<double> delta_uI(N_PARAM);
			//double delta_phiR = 0;
			//double delta_phiI = 0;

			//cout << "solve JR" << endl;
			//WriteDataToFile(JR, "BBJR_" + to_string(iterations), "JR");
			//WriteDataToFile(gR_n, "BBgR_n_" + to_string(iterations), "gR_n");
			//WriteDataToFile(uR_n, "BBuR_n_" + to_string(iterations), "uR_n");
			SolveEquationSystemArmadillo(JR, gR_n, delta_uR);
			//WriteDataToFile(delta_uR, "BBdelta_uR_" + to_string(iterations), "delta_uR");

			//cout << "solve JI" << endl;
			//WriteDataToFile(JI, "BBJI_" + to_string(iterations), "JI");
			//WriteDataToFile(gI_n, "BBgI_n_" + to_string(iterations), "gI_n");
			//WriteDataToFile(uI_n, "BBuI_n_" + to_string(iterations), "uI_n");
			SolveEquationSystemArmadillo(JI, gI_n, delta_uI);
			//WriteDataToFile(delta_uI, "BBdelta_uI_" + to_string(iterations), "delta_uI");

			double max = abs(delta_uR[0] / uR_n[0]);
			double maxNoAbs = delta_uR[0] / uR_n[0];
			int maxIndex = 0;
			double tmp = 0;
			for (int i = 0; i < N_PARAM; i++)
			{
				tmp = abs(delta_uR[i] / uR_n[i]);
				if (tmp > max)
				{
					max = tmp;
					maxNoAbs = delta_uR[i] / uR_n[i];
					maxIndex = i;
				}
				tmp = abs(delta_uI[i] / uI_n[i]);
				if (tmp > max)
				{
					maxNoAbs = delta_uI[i] / uI_n[i];
					max = tmp;
					maxIndex = -i;
				}
			}
			diff = max;
			cout << "diff=" << diff << " (" << maxNoAbs << "@" << maxIndex << ")" << endl;

			cout << "compute u_" << iterations << endl;
			uR_n = uR_n - delta_uR;
			uI_n = uI_n - delta_uI;
		}

		MPIMethods::BroadcastValue(&diff);
		BroadcastNewParameters(uR_n, uI_n, &phiR_n, &phiI_n);
	}

	if (processRank == rootRank)
	{
		//TODO: check calculation of phiR and phII (but should not affect results since phiR is set in NormalizeWavefunction() and phiI is only global phase factor...)
		vector<double> final_uDotR;
		vector<double> final_uDotI;
		double final_phiDotR;
		double final_phiDotI;
		final_uDotR = (uR - uR_n) / dt;
		final_uDotI = (uI - uI_n) / dt;
		CalculatePhiDot(final_uDotR, final_uDotI, &final_phiDotR, &final_phiDotI);

		uR = uR_n;
		uI = uI_n;
		*phiR = *phiR + final_phiDotR * dt;
		*phiI = *phiI + final_phiDotI * dt;
	}
}

void CalculateNextParametersCrankNicolson(double dt, vector<double>& uR, vector<double>& uI, double *phiR, double *phiI)
{

	vector<double> uR_Start(uR);
	vector<double> uI_Start(uI);
	double phiR_Start = *phiR;
	double phiI_Start = *phiI;

	vector<double> uDotR(N_PARAM);
	vector<double> uDotI(N_PARAM);
	double phiDotR = 0;
	double phiDotI = 0;

	vector<double> uDotR_n(N_PARAM);
	vector<double> uDotI_n(N_PARAM);
	double phiDotR_n = 0;
	double phiDotI_n = 0;

	vector<double> uDotR_n_delta(N_PARAM);
	vector<double> uDotI_n_delta(N_PARAM);
	double phiDotR_n_delta = 0;
	double phiDotI_n_delta = 0;

	vector<vector<double> > JR(N_PARAM);
	vector<vector<double> > JI(N_PARAM);
	for (int i = 0; i < N_PARAM; i++)
	{
		JR[i].resize(N_PARAM);
		JI[i].resize(N_PARAM);
	}

	if (processRank == rootRank)
	{
		SolveForParametersDot(uDotR, uDotI, &phiDotR, &phiDotI);

		uR_Start = uR_Start + uDotR * dt;
		uI_Start = uI_Start + uDotI * dt;
		phiR_Start = phiR_Start + phiDotR * dt;
		phiI_Start = phiI_Start + phiDotI * dt;
	}

	BroadcastNewParameters(uR_Start, uI_Start, &phiR_Start, &phiI_Start);
	BroadcastNewParameters(uDotR, uDotI, &phiDotR, &phiDotI);

	double eps = 1e-2;
	double diff = 1;
	double delta = 1e-4;
	int iterations = 0;

	vector<double> uR_n(uR_Start);
	vector<double> uI_n(uI_Start);
	double phiR_n = phiR_Start;
	double phiI_n = phiI_Start;
	//vector<double> uR_n(uR);
	//vector<double> uI_n(uI);
	//double phiR_n = *phiR;
	//double phiI_n = *phiI;

	//vector<double> directionR(N_PARAM);
	//vector<double> directionI(N_PARAM);
	//for (int i = 0; i < N_PARAM; i++)
	//{
	//	directionR[i] = uDotR[i] > 0 ? 1 : -1;
	//	directionI[i] = uDotI[i] > 0 ? 1 : -1;
	//}

	while (diff > eps && iterations < 20)
	{
		iterations++;
		ParallelUpdateExpectationValuesForGivenSamples(correlatedSamplingData, uR_n, uI_n, phiR_n, phiI_n, true);
		if (processRank == rootRank)
		{
			SolveForParametersDot(uDotR_n, uDotI_n, &phiDotR_n, &phiDotI_n);
		}
		BroadcastNewParameters(uDotR_n, uDotI_n, &phiDotR_n, &phiDotI_n);

		vector<double> gR_n = uR_n - uR - uDotR_n * dt; //TODO: gR_n and uDotR_n is only used by rootRank?!?
		vector<double> gI_n = uI_n - uI - uDotI_n * dt;

		for (int i = 0; i < N_PARAM; i++)
		{
			double dUR = delta * uR[i]; // * directionR[i];
			uR_n[i] += dUR;
			ParallelUpdateExpectationValuesForGivenSamples(correlatedSamplingData, uR_n, uI_n, phiR_n, phiI_n, true);
			if (processRank == rootRank)
			{
				SolveForParametersDot(uDotR_n_delta, uDotI_n_delta, &phiDotR_n_delta, &phiDotI_n_delta);
				JR[i] = (uDotR_n_delta - uDotR_n) * dt / dUR * (-1.0);
				JR[i][i] = JR[i][i] + 1.0;
			}
			uR_n[i] -= dUR;

			double dUI = delta * uI[i]; // * directionI[i];
			uI_n[i] += dUI;
			ParallelUpdateExpectationValuesForGivenSamples(correlatedSamplingData, uR_n, uI_n, phiR_n, phiI_n, true);
			if (processRank == rootRank)
			{
				SolveForParametersDot(uDotR_n_delta, uDotI_n_delta, &phiDotR_n_delta, &phiDotI_n_delta);
				JI[i] = (uDotI_n_delta - uDotI_n) * dt / dUI * (-1.0);
				JI[i][i] = JI[i][i] + 1.0;
			}
			uI_n[i] -= dUI;
		}
		if (processRank == rootRank)
		{
			vector<double> delta_uR(N_PARAM);
			vector<double> delta_uI(N_PARAM);
			//double delta_phiR = 0;
			//double delta_phiI = 0;

			cout << "solve JR" << endl;
			//WriteDataToFile(JR, "BBJR_" + to_string(iterations), "JR");
			//WriteDataToFile(gR_n, "BBgR_n_" + to_string(iterations), "gR_n");
			//WriteDataToFile(uR_n, "BBuR_n_" + to_string(iterations), "uR_n");
			SolveEquationSystemArmadillo(JR, gR_n, delta_uR);
			//WriteDataToFile(delta_uR, "BBdelta_uR_" + to_string(iterations), "delta_uR");

			cout << "solve JI" << endl;
			//WriteDataToFile(JI, "BBJI_" + to_string(iterations), "JI");
			//WriteDataToFile(gI_n, "BBgI_n_" + to_string(iterations), "gI_n");
			//WriteDataToFile(uI_n, "BBuI_n_" + to_string(iterations), "uI_n");
			SolveEquationSystemArmadillo(JI, gI_n, delta_uI);
			//WriteDataToFile(delta_uI, "BBdelta_uI_" + to_string(iterations), "delta_uI");

			double max = abs(delta_uR[0] / uR_n[0]);
			double maxNoAbs = delta_uR[0] / uR_n[0];
			int maxIndex = 0;
			double tmp = 0;
			for (int i = 0; i < N_PARAM; i++)
			{
				tmp = abs(delta_uR[i] / uR_n[i]);
				if (tmp > max)
				{
					max = tmp;
					maxNoAbs = delta_uR[i] / uR_n[i];
					maxIndex = i;
				}
				tmp = abs(delta_uI[i] / uI_n[i]);
				if (tmp > max)
				{
					maxNoAbs = delta_uI[i] / uI_n[i];
					max = tmp;
					maxIndex = -i;
				}
			}
			diff = max;
			cout << "diff=" << diff << " (" << maxNoAbs << "@" << maxIndex << ")" << endl;

			cout << "compute u_n+1" << endl;
			uR_n = uR_n - delta_uR;
			uI_n = uI_n - delta_uI;
		}

		MPIMethods::BroadcastValue(&diff);
		BroadcastNewParameters(uR_n, uI_n, &phiR_n, &phiI_n);
	}

	if (processRank == rootRank)
	{
		//TODO: check calculation of phiR and phII (but should not affect results since phiR is set in NormalizeWavefunction() and phiI is only global phase factor...)
		vector<double> final_uDotR;
		vector<double> final_uDotI;
		double final_phiDotR;
		double final_phiDotI;
		final_uDotR = ((uR - uR_n) / dt + uDotR) / 2.0; //INFO: mean value of backward euler and forward euler
		final_uDotI = ((uI - uI_n) / dt + uDotI) / 2.0;
		CalculatePhiDot(final_uDotR, final_uDotI, &final_phiDotR, &final_phiDotI);

		uR = uR + final_uDotR * dt;
		uI = uI + final_uDotI * dt;
		*phiR = *phiR + final_phiDotR * dt;
		*phiI = *phiI + final_phiDotI * dt;
	}
}

void CalculateNextParameters(double dt, vector<vector<double> >& R, vector<double>& uR, vector<double>& uI, double *phiR, double *phiI)
{
	Timer t;
	if (processRank == rootRank)
	{
		t.start();
	}

	if (ODE_SOLVER_TYPE == 0)
	{
		CalculateNextParametersEuler(dt, uR, uI, phiR, phiI);
	}
	else if (ODE_SOLVER_TYPE == 1)
	{
		vector<vector<double> >& Rcopy(R);
		CalculateNextParametersPC(dt, Rcopy, uR, uI, phiR, phiI);
	}
	else if (ODE_SOLVER_TYPE == 2)
	{
		vector<vector<double> >& Rcopy(R);
		CalculateNextParametersRK4(dt, Rcopy, uR, uI, phiR, phiI);
	}
	else if (ODE_SOLVER_TYPE == 3)
	{
		CalculateNextParametersRK4ReuseSamples(dt, uR, uI, phiR, phiI);
	}
	else if (ODE_SOLVER_TYPE == 4)
	{
		if (sys->GetStep() < 2)
		{
			CalculateNextParametersRK4ReuseSamples(dt, uR, uI, phiR, phiI);
		}
		else
		{
			CalculateNextParametersImplicitEuler(dt, uR, uI, phiR, phiI);
		}
	}
	else if (ODE_SOLVER_TYPE == 5)
	{
		CalculateNextParametersPCReuseSamples(dt, uR, uI, phiR, phiI);
	}
	else if (ODE_SOLVER_TYPE == 6)
	{
		if (sys->GetStep() < 2)
		{
			CalculateNextParametersRK4ReuseSamples(dt, uR, uI, phiR, phiI);
		}
		else
		{
			CalculateNextParametersCrankNicolson(dt, uR, uI, phiR, phiI);
		}
	}

	if (processRank == rootRank)
	{
		t.stop();
		Log("DGL duration = " + to_string(t.duration()) + " ms");
	}
}

void CalculateNextParameters(double dt, vector<vector<double> >& R, vector<double>& uR, vector<double>& uI)
{
	double tmpPhiR = 0;
	double tmpPhiI = 0;
	CalculateNextParameters(dt, R, uR, uI, &tmpPhiR, &tmpPhiI);
}

////////////
/// misc ///
////////////

void NormalizeWavefunction(double wf, double *phiR)
{
	*phiR = *phiR - log(wf);
}

void AdjustParameters(vector<double>& uR, vector<double>& uI, double *phiR, double *phiI)
{
	//double last = uR.back() + 0.00908425;
	//for (unsigned int i = 0; i < uR.size(); i++)
	//{
	//	uR[i] += -last;
	//}
	//uR[uR.size() - 1] = 0.0;
	//uR[0] = uR[2];
	//uR[1] = uR[2];

	//for (int i = 40; i < N_PARAM; i++)
	//{
	//	//uR[i] *= (1.0 / (1.0 + exp((i-43.0)/1.3)) + 3.0) / 4.0;
	//	//uR[i] = (uR[i] + (i - 35.0) * uRList.back()[i]) / (1.0 + i - 35.0);
	//	//uR[i] += -uR[45];
	//}

	//for (int i = 91; i < N_PARAM; i++)
	//{
	//	uR[i] = -uR[90];
	//}

	//for (int i = N_PARAM - 15; i < N_PARAM; i++)
	//{
	//	uR[i] = (uR[i] + (i - 14) * uRList.back()[i]) / (1.0 + (i - 14));
	//	uI[i] = (uI[i] + (i - 14) * uIList.back()[i]) / (1.0 + (i - 14));
	//}

	for (int i = 0; i < N_PARAM; i++)
	{
		//uR[i] = (uR[i] + 0.1 * i * uRList.back()[i]) / (1.0 + 0.1 * i);
		//uI[i] = (uI[i] + 0.1 * i * uIList.back()[i]) / (1.0 + 0.1 * i);
		uR[i] = uR[i] - uR[N_PARAM - 1];
		uI[i] = uI[i] - uI[N_PARAM - 1];
	}
}

void AlignCoordinates(vector<vector<double> >& R)
{
	if (sys->USE_NIC)
	{
		MoveCoordinatesToFirstCell(R);
	}
	else
	{
		if (sys->USE_MOVE_COM_TO_ZERO)
		{
			MoveCenterOfMassToZero(R);
		}
	}
}

//TODO: find better criteria when to accept new parameters
bool AcceptNewParams(vector<double>& uR, vector<double>& uI, double phiR, double phiI, double timestep)
{
	bool accept = true;
	if (doNotAcceptStep) //INFO: An error already occured in PerformCholeskyDecomposition
	{
		doNotAcceptStep = false;
		return false;
	}
	else if (PARAMETER_ACCEPTANCE_CHECK_TYPE == 0)
	{
		return true;
	}
	else if (PARAMETER_ACCEPTANCE_CHECK_TYPE == 1)
	{
		if (uRList.size() > 1)
		{
			for (int p = 0; p < N_PARAM; p++)
			{
				if (!isfinite(uR[p]))
				{
					return false;
				}
			}
		}
	}
	else if (PARAMETER_ACCEPTANCE_CHECK_TYPE == 2)
	{
		double threshold = 0.1;
		if (uRList.size() > 1)
		{
			for (int p = 0; p < N_PARAM; p++)
			{
				if (!isfinite(uR[p]))
				{
					return false;
				}
				if (abs(1.0 - uR[p] / uRList.back()[p]) > threshold)
				{
					return false;
				}
			}
		}
	}
	else if (PARAMETER_ACCEPTANCE_CHECK_TYPE == 3)
	{
		unsigned int lastNValues = 10;
		double value = 0;
		double newValue = 0;
		if (uRList.size() > lastNValues + 1)
		{
			for (int p = 0; p < N_PARAM; p++)
			{
				if (!isfinite(uR[p]))
				{
					return false;
				}
				value = 0;
				for (unsigned int i = 0; i < lastNValues; i++)
				{
					value += abs(uRListDiffs[uRListDiffs.size() - i - 1][p]);
				}
				newValue = abs(uR[p] - uRList.back()[p]);
				if (newValue > value * 5)
				{
					return false;
				}
			}
		}
	}
	else if (PARAMETER_ACCEPTANCE_CHECK_TYPE == 4)
	{
		double thresholdParamR = timestep * 1e5;
		double thresholdEnergyR = timestep * 1e5;
		//double thresholdEnergyI = 1.0;
		if (AllLocalEnergyR.size() > 1 && abs(GetRelativeErrorLastElements(AllLocalEnergyR)) > thresholdEnergyR)
		{
			cout << "Energy change: " << GetRelativeErrorLastElements(AllLocalEnergyR) << endl;
			return false;
		}
		//if (AllLocalEnergyI.size() > 5 && abs(GetRelativeErrorLastElements(AllLocalEnergyI)) > thresholdEnergyI)
		//{
		//	PrintData(AllLocalEnergyI);
		//	cout << "err=" << GetRelativeErrorLastElements(AllLocalEnergyI) << endl;
		//	return false;
		//}
		if (uRList.size() > 1)
		{
			for (int p = 0; p < N_PARAM; p++)
			{
				if (!isfinite(uR[p]))
				{
					return false;
				}
				if (abs(GetRelativeError(uRList.back()[p], uR[p])) > thresholdParamR)
				{
					cout << "Parameter nr " << p << " change: " << GetRelativeError(uRList.back()[p], uR[p]) << endl;
					return false;
				}
			}
		}
	}
	else if (PARAMETER_ACCEPTANCE_CHECK_TYPE == 5)
	{
		if (uRList.size() > 1)
		{
			for (int p = 0; p < N_PARAM; p++)
			{
				if (!isfinite(uR[p]))
				{
					return false;
				}
			}
		}
		if (AllLocalEnergyR.size() > 20)
		{
			double maxLastRelativeDiff = 0.0;
			double relativeDiff;
			for (int i = 0; i < 10; i++)
			{
				maxLastRelativeDiff = max(maxLastRelativeDiff, abs(GetRelativeError(AllLocalEnergyR[AllLocalEnergyR.size() - i - 3], AllLocalEnergyR[AllLocalEnergyR.size() - 2])));
			}
			relativeDiff = abs(GetRelativeError(AllLocalEnergyR[AllLocalEnergyR.size() - 2], AllLocalEnergyR.back()));
			if (relativeDiff / 10.0 > maxLastRelativeDiff)
			{
				return false;
			}
		}
	}
	return accept;
}

void PerformSystemParameterChange()
{
	string errors;
	bool successful;
	Json::Value configData;
	Json::CharReaderBuilder builder;
	ifstream configFile("./param", ifstream::binary);

	successful = Json::parseFromStream(builder, configFile, &configData, &errors);

	if (successful)
	{
		for (auto ci : configItems)
		{
			if (ci.allowChangeAtRuntime && !configData[ci.name].empty())
			{
				string oldValue = ci.getJsonString();
				ci.setValue(configData[ci.name]);
				Log("changed: " + ci.name + " from \"" + oldValue + "\" to \"" + ci.getJsonString() + "\"", WARNING);
			}
		}
	}
	else
	{
		Log("could not read param change file", ERROR);
	}
}

void SetBackNSteps(int n, vector<vector<double> >& R, vector<double>& uR, vector<double>& uI, double phiR, double phiI, int& step, double dynTimeStep)
{
	//TODO: also check for urListDiff, correlated sampling, ...
	if (processRank == rootRank)
	{
		for (int i = 0; i < n; i++)
		{
			currentTime -= dynTimeStep;
			step--;
			sys->SetTime(currentTime);
			sys->SetStep(step);
			times.pop_back();

			AllLocalEnergyR.pop_back();
			AllLocalEnergyI.pop_back();
			AllLocalOperators.pop_back();
			AllOtherExpectationValues.pop_back();
			AllParametersR.pop_back();
			AllParametersI.pop_back();
			uR = uRList.back();
			uI = uIList.back();
			phiR = phiRList.back();
			phiI = phiIList.back();
		}
	}
	AlignCoordinates(R);
	BroadcastNewParameters(uR, uI, &phiR, &phiI);
	sys->CalculateWavefunction(R, uR, uI, phiR, phiI);
	for (int i = 0; i < MC_VERY_FIRST_NINITIALIZATIONSTEPS; i++)
	{
		DoMetropolisStep(R, uR, uI, phiR, phiI);
	}
	sys->CalculateWavefunction(R, uR, uI, phiR, phiI);
}

bool InitializePhysicalSystem()
{
	bool successful = true;
	if (SYSTEM_TYPE == "HeDrop")
	{
		sys = new PhysicalSystems::HeDrop(SYSTEM_PARAMS, configDirectory);
	}
	else if (SYSTEM_TYPE == "HeBulk")
	{
		sys = new PhysicalSystems::HeBulk(SYSTEM_PARAMS, configDirectory);
	}
	else if (SYSTEM_TYPE == "BulkSplines")
	{
		//sys = new BulkSplines(SYSTEM_PARAMS, configDirectory);
		sys = new PhysicalSystems::BulkSplinesPhi(SYSTEM_PARAMS, configDirectory);
	}
	else if (SYSTEM_TYPE == "BulkSplinesScaled")
	{
		sys = new PhysicalSystems::BulkSplinesScaled(SYSTEM_PARAMS, configDirectory);
	}
	else if (SYSTEM_TYPE == "BulkQT")
	{
		//sys = new BulkQT(SYSTEM_PARAMS, configDirectory);
		sys = new PhysicalSystems::BulkQTPhi(SYSTEM_PARAMS, configDirectory);
	}
	else if (SYSTEM_TYPE == "GaussianWavepacket")
	{
		sys = new PhysicalSystems::GaussianWavepacket(SYSTEM_PARAMS, configDirectory);
	}
	else if (SYSTEM_TYPE == "BosonsBulk")
	{
		sys = new PhysicalSystems::BosonsBulk(SYSTEM_PARAMS, configDirectory);
	}
	else if (SYSTEM_TYPE == "BosonsBulkDamped")
	{
		sys = new PhysicalSystems::BosonsBulkDamped(SYSTEM_PARAMS, configDirectory);
	}
	else if (SYSTEM_TYPE == "BosonCluster")
	{
		sys = new PhysicalSystems::BosonCluster(SYSTEM_PARAMS, configDirectory);
		if (USE_NURBS == 1)
		{
			dynamic_cast<PhysicalSystems::BosonCluster*>(sys)->SetNodes(NURBS_GRID);
		}
		dynamic_cast<PhysicalSystems::BosonCluster*>(sys)->SetDensityProfileBinCount(GR_BIN_COUNT);
	}
	else if (SYSTEM_TYPE == "BosonClusterWithLog")
	{
		sys = new PhysicalSystems::BosonClusterWithLog(SYSTEM_PARAMS, configDirectory);
		if (USE_NURBS == 1)
		{
			dynamic_cast<PhysicalSystems::BosonClusterWithLog*>(sys)->SetNodes(NURBS_GRID);
		}
		dynamic_cast<PhysicalSystems::BosonClusterWithLog*>(sys)->SetDensityProfileBinCount(GR_BIN_COUNT);
	}
	else if (SYSTEM_TYPE == "BosonClusterWithLogParam")
	{
		sys = new PhysicalSystems::BosonClusterWithLogParam(SYSTEM_PARAMS, configDirectory);
		if (USE_NURBS == 1)
		{
			dynamic_cast<PhysicalSystems::BosonClusterWithLogParam*>(sys)->SetNodes(NURBS_GRID);
		}
		dynamic_cast<PhysicalSystems::BosonClusterWithLogParam*>(sys)->SetDensityProfileBinCount(GR_BIN_COUNT);
	}
	else if (SYSTEM_TYPE == "BosonMixtureCluster")
	{
		sys = new PhysicalSystems::BosonMixtureCluster(SYSTEM_PARAMS, configDirectory);
		if (USE_NURBS == 1)
		{
			dynamic_cast<PhysicalSystems::BosonMixtureCluster*>(sys)->SetNodes(NURBS_GRID);
		}
		dynamic_cast<PhysicalSystems::BosonMixtureCluster*>(sys)->SetDensityProfileBinCount(GR_BIN_COUNT);
		dynamic_cast<PhysicalSystems::BosonMixtureCluster*>(sys)->SetParticleType(PARTICLE_TYPES);
	}
	else if (SYSTEM_TYPE == "HardSphereBosons")
	{
		sys = new PhysicalSystems::HardSphereBosons(SYSTEM_PARAMS, configDirectory);
	}
	else if (SYSTEM_TYPE == "HardSphereBosonsExp")
	{
		sys = new PhysicalSystems::HardSphereBosonsExp(SYSTEM_PARAMS, configDirectory);
	}
	else if (SYSTEM_TYPE == "NUBosonsBulk")
	{
		sys = new PhysicalSystems::NUBosonsBulk(SYSTEM_PARAMS, configDirectory);
		if (USE_NURBS == 1)
		{
			dynamic_cast<PhysicalSystems::NUBosonsBulk*>(sys)->SetNodes(NURBS_GRID);
		}
		dynamic_cast<PhysicalSystems::NUBosonsBulk*>(sys)->SetGrBinCount(GR_BIN_COUNT);
	}
	else if (SYSTEM_TYPE == "PBosonsBulk")
	{
		sys = new PhysicalSystems::PBosonsBulk(SYSTEM_PARAMS, configDirectory);
		if (USE_NURBS == 1)
		{
			dynamic_cast<PhysicalSystems::PBosonsBulk*>(sys)->SetNodes(NURBS_GRID);
		}
	}
	else
	{
		if (processRank == rootRank)
		{
			Log("System type \"" + SYSTEM_TYPE + "\" not available.", ERROR);
		}
		successful = false;
	}
	return successful;
}

////////////
/// main ///
////////////

int mainMPI(int argc, char** argv)
{
	//////////////////////////
	/// Init MPI processes ///
	//////////////////////////

	char processName[80];
	int processNameLength;

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &processRank);
	MPI_Comm_size(MPI_COMM_WORLD, &numOfProcesses);
	MPI_Get_processor_name(processName, &processNameLength);

	MPIMethods::numOfProcesses = numOfProcesses;
	MPIMethods::processRank = processRank;
	MPIMethods::rootRank = rootRank;

	MPIMethods::GetCPUAllocation(true);

	////////////////////////////////////////
	/// Read and Broadcast configuration ///
	////////////////////////////////////////

	RegisterAllConfigItems();
	if (processRank == rootRank)
	{
		Log("Master process started...");
		cout << "read config from path=\"" << configFilePath << "\"" << endl;

		if (FileExist(configFilePath))
		{
			ReadConfig(configFilePath);
			if (requiredConfigVersion != CONFIG_VERSION)
			{
				Log("Config file is out of date.\nrequired: \"" + requiredConfigVersion + "\"\nfound: \"" + CONFIG_VERSION + "\"", ERROR);
				return 1;
			}
			CreateOutputDirectory();
			//PrintConfig();
		}
		else
		{
			cout << "config file not found. path=\"" << configFilePath << "\"" << endl;
			MPI_Abort(MPI_COMM_WORLD, 0);
			return 1;
		}
	}
	configDirectory = configFilePath.substr(0, configFilePath.find_last_of("\\/")) + "/"; //TODO: restrict all file access to main process
	MPIMethods::BroadcastValue(&OUT_DIR, 300);
	//cout << "configDirectory=\"" << configDirectory << "\"" << endl;
	BroadcastConfig();
	//PrintConfig();

	//////////////////////////////////
	/// Initialize Physical System ///
	//////////////////////////////////

	if (!InitializePhysicalSystem())
	{
		return 1;
	}

	/////////////////////////////
	/// other Initializations ///
	/////////////////////////////

	Init();
	if (processRank == rootRank)
	{
		//cout << "uniform: " << random01() << endl << "normal: " << randomNormal() << endl;
		//cout << "uniform: " << random01() << endl << "normal: " << randomNormal() << endl;
		//cout << "uniform: " << random01() << endl << "normal: " << randomNormal() << endl;
	}

	vector<vector<double> > R(N);
	vector<double> uR(N_PARAM);
	vector<double> uI(N_PARAM);
	double phiR;
	double phiI;
	PARAMS_REAL.resize(N_PARAM, 0.0);
	PARAMS_IMAGINARY.resize(N_PARAM, 0.0);

	for (unsigned int i = 0; i < R.size(); i++)
	{
		R[i].resize(DIM);
	}

	InitCoordinateConfiguration(R);
	if (processRank == rootRank)
	{
		WriteParticleInputFile("particleconfiguration", R);
	}
	if (processRank == rootRank)
	{
		for (int i = 0; i < N_PARAM; i++)
		{
			uR[i] = PARAMS_REAL[i];
			uI[i] = PARAMS_IMAGINARY[i];
		}
		phiR = PARAM_PHIR;
		phiI = PARAM_PHII;
		//PrintData(uR);
		//PrintData(uI);

		if (processRank == rootRank)		// && USE_ADJUST_PARAMETERS == true)
		{
			//AdjustParameters(uR, uI, &phiR, &phiI);
		}

		WriteDataToFile(uR, "parametersR0", "parameterR");
		WriteDataToFile(uI, "parametersI0", "parameterI");

		WriteDataToFile(AllLocalEnergyR, "AllLocalEnergyR", "ER", WRITE_EVERY_NTH_STEP_TO_FILE);
		WriteDataToFile(AllLocalEnergyI, "AllLocalEnergyI", "EI", WRITE_EVERY_NTH_STEP_TO_FILE);
		WriteDataToFile(AllLocalOperators, "AllLocalOperators", "<O_k>", WRITE_EVERY_NTH_STEP_TO_FILE);
		WriteDataToFile(AllOtherExpectationValues, "AllOtherExpectationValues", "kinetic, potential, wf, g(r)", WRITE_EVERY_NTH_STEP_TO_FILE);
		WriteDataToFile(AllParametersR, "AllParametersR", "uR", WRITE_EVERY_NTH_STEP_TO_FILE);
		WriteDataToFile(AllParametersI, "AllParametersI", "uI", WRITE_EVERY_NTH_STEP_TO_FILE);
	}

	sys->InitSystem();
	PostSystemInit();

	////////////////////////
	/// start simulation ///
	////////////////////////

	Timer t;
	int step = 0;
	int acceptNewParams;
	int nrOfAcceptParameterTrials;
	int maxNrOfAcceptParameterTrials = 2;

	double nrOfSamplesToUpdate;
	double percentForExtraSampleToUpdate;
	percentForExtraSampleToUpdate = modf(MC_NSTEPS * UPDATE_SAMPLES_PERCENT / 100.0, &nrOfSamplesToUpdate);

	double dynTimestep = TIMESTEP;
	AlignCoordinates(R);
	BroadcastNewParameters(uR, uI, &phiR, &phiI);
	sys->CalculateWavefunction(R, uR, uI, phiR, phiI);
	for (int i = 0; i < MC_VERY_FIRST_NINITIALIZATIONSTEPS; i++)
	{
		DoMetropolisStep(R, uR, uI, phiR, phiI);
	}
	sys->CalculateWavefunction(R, uR, uI, phiR, phiI);
	for (currentTime = 0; currentTime <= TOTALTIME; currentTime += dynTimestep)
	{
		if (processRank == rootRank)
		{
			t.start();
		}
		acceptNewParams = 0;
		nrOfAcceptParameterTrials = 0;
		MC_NSTEPS = mc_nsteps_original;
		nAcceptances = 0;
		nTrials = 0;
		step++;
		sys->SetTime(currentTime);
		sys->SetStep(step);
		times.push_back(currentTime);

		BroadcastNewParameters(uR, uI, &phiR, &phiI);
		if (processRank == rootRank)
		{
			if (uRList.size() > 0)
			{
				uRListDiffs.push_back(uRList.back() - uR);
				uIListDiffs.push_back(uIList.back() - uI);
				phiRListDiffs.push_back(phiRList.back() - phiR);
				phiIListDiffs.push_back(phiIList.back() - phiI);
			}
			uRList.push_back(uR);
			uIList.push_back(uI);
			phiRList.push_back(phiR);
			phiIList.push_back(phiI);
			if (uRList.size() > 105) //INFO: always keep the last few parameters in a list. at the end of the simulation the average of the last steps is calculated for a final value of the parameter
			{
				uRList.erase(uRList.begin());
				uIList.erase(uIList.begin());
				phiRList.erase(phiRList.begin());
				phiIList.erase(phiIList.begin());
				uRListDiffs.erase(uRListDiffs.begin());
				uIListDiffs.erase(uIListDiffs.begin());
				phiRListDiffs.erase(phiRListDiffs.begin());
				phiIListDiffs.erase(phiIListDiffs.begin());
			}
		}
		//PrintParameters(uR);
		//PrintParameters(uI);
		//cout << "phiR=" << phiR << ", phiI=" << phiI << endl;

		////////////////////////////////////////
		/// try to calculate next parameters ///
		////////////////////////////////////////

		while (acceptNewParams != 1 && nrOfAcceptParameterTrials < maxNrOfAcceptParameterTrials)
		{
			nrOfAcceptParameterTrials++;
			MC_NSTEPS *= nrOfAcceptParameterTrials;
			AlignCoordinates(R);

			if (step < 2 || UPDATE_SAMPLES_EVERY_NTH_STEP == 0)
			{
				ParallelUpdateExpectationValues(R, uR, uI, phiR, phiI);
				//if (ODE_SOLVER_TYPE == 4 || ODE_SOLVER_TYPE == 6)
				//{
				//	dynTimestep = 0.01 * TIMESTEP;
				//}
			}
			else
			{
				dynTimestep = TIMESTEP;
				if ((UPDATE_SAMPLES_EVERY_NTH_STEP > 0 && step % UPDATE_SAMPLES_EVERY_NTH_STEP == 0) || nrOfAcceptParameterTrials > 0)
				{
					//UpdateSamplesRandom(nrOfSamplesToUpdate, correlatedSamplingData, uR, uI, phiR, phiI);
					////cout << "nrOfSamplesToUpdate " << nrOfSamplesToUpdate << endl;
					////cout << "percentForExtraSampleToUpdate " << percentForExtraSampleToUpdate << endl;
					//if (random01() < percentForExtraSampleToUpdate)
					//{
					//	//cout << "update extra sample..." << endl;
					//	UpdateSamplesRandom(1, correlatedSamplingData, uR, uI, phiR, phiI);
					//}
					UpdateSamplesConsecutive(nrOfSamplesToUpdate, correlatedSamplingData, uR, uI, phiR, phiI);
				}
				ParallelUpdateExpectationValuesForGivenSamples(correlatedSamplingData, uR, uI, phiR, phiI);

				//bool update = NeedToUpdateSamples(mcSamples, uR, uI, phiR, phiI);
				//vector<bool> allUpdateValues = MPIMethods::CollectValues(update);
				//if (processRank == rootRank)
				//{
				//	PrintData(allUpdateValues);
				//}
				//if (MPIMethods::IsAnyFalse(update))
				//{
				//	if (update)
				//	{
				//		UpdateAllSamples(mcSamples, uR, uI, phiR, phiI);
				//	}
				//	ParallelUpdateExpectationValuesForGivenSamples(mcSamples, uR, uI, phiR, phiI);
				//	dynTimestep = TIMESTEP;// * 10.0;
				//}
				//else
				//{
				//	ParallelUpdateExpectationValues(R, uR, uI, phiR, phiI);
				//	dynTimestep = TIMESTEP;
				//}

				//if (MPIMethods::IsAnyTrue(update))
				//{
				//	ParallelUpdateExpectationValues(R, uR, uI, phiR, phiI);
				//}
				//else
				//{
				//	if (processRank == rootRank)
				//	{
				//		Log("reuse samples!", WARNING);
				//	}
				//	ParallelUpdateExpectationValuesForGivenSamples(mcSamples, uR, uI, phiR, phiI);
				//}
			}
			//if (step < 2)
			//{
			//	ParallelUpdateExpectationValues(R, uR, uI, phiR, phiI);
			//}
			//else
			//{
			//	//if (step % 2 == 1)
			//	{
			//		UpdateSamples(mcSamples, uR, uI, phiR, phiI);
			//	}
			//	ParallelUpdateExpectationValuesForGivenSamples(mcSamples, uR, uI, phiR, phiI);
			//}
			//ParallelUpdateExpectationValues(R, uR, uI, phiR, phiI);

			double avgAcceptances = MPIMethods::ReduceToAverage(&nAcceptances);
			if (processRank == rootRank)
			{
				if (WRITE_SINGLE_FILES == 1 && step % WRITE_EVERY_NTH_STEP_TO_FILE == 0)
				{
					WriteDataToFile(localOperators, "localOperators" + to_string(step), "localOperators");
					WriteDataToFile(localEnergyR, "localEnergyR" + to_string(step), "localEnergyR");
					WriteDataToFile(localEnergyI, "localEnergyI" + to_string(step), "localEnergyI");
					WriteDataToFile(localOperatorsMatrix, "localOperatorsMatrix" + to_string(step), "localOperatorsMatrix");
					WriteDataToFile(localOperatorlocalEnergyR, "localOperatorlocalEnergyR" + to_string(step), "localOperatorlocalEnergyR");
					WriteDataToFile(localOperatorlocalEnergyI, "localOperatorlocalEnergyI" + to_string(step), "localOperatorlocalEnergyI");
					WriteDataToFile(otherExpectationValues, "otherExpectationValues" + to_string(step), "Ekin, Ekin_cor, Epot, Epot_corr, wf, g(r)_1, ..., g(r)_100");
				}
				cout << "t=" << currentTime << endl;
				cout << "Acceptance AVG: " << (avgAcceptances / (nTrials / 100.0)) << "% (" << avgAcceptances << "/" << nTrials << ")" << endl;
				//cout << "localEnergyR=" << localEnergyR << " (" << otherExpectationValues[0] << " + " << otherExpectationValues[1] << ")" << endl;
				cout << "localEnergyR/N=" << localEnergyR / (double) N << " (" << otherExpectationValues[0] / (double) N << " + " << otherExpectationValues[1] / (double) N << " + " <<
				//otherExpectationValues[2] / (double) N << " + " <<
				//otherExpectationValues[3] / (double) N << ")" <<
						endl;
				//cout << "localEnergyR/N=" << localEnergyR / (double)N << endl;
				AllLocalEnergyR.push_back(localEnergyR);
				AllLocalEnergyI.push_back(localEnergyI);
				AllLocalOperators.push_back(localOperators);
				AllOtherExpectationValues.push_back(otherExpectationValues);
				AllParametersR.push_back(uR);
				AllParametersR[AllParametersR.size() - 1].push_back(phiR);
				AllParametersI.push_back(uI);
				AllParametersI[AllParametersI.size() - 1].push_back(phiI);
			}
			if (sys->USE_NORMALIZATION_AND_PHASE)
			{
				CalculateNextParameters(dynTimestep, R, uR, uI, &phiR, &phiI);
				if (processRank == rootRank && USE_NORMALIZE_WF == 1)
				{
					NormalizeWavefunction(sys->GetWf(), &phiR);
					//or better use
					//NormalizeWavefunction(otherExpectationValues[2], &phiR);
				}
			}
			else
			{
				CalculateNextParameters(dynTimestep, R, uR, uI);
			}
			if (processRank == rootRank && USE_ADJUST_PARAMETERS == 1)
			{
				AdjustParameters(uR, uI, &phiR, &phiI);
			}
			if (USE_PARAMETER_ACCEPTANCE_CHECK == 1)
			{
				if (processRank == rootRank)
				{
					acceptNewParams = AcceptNewParams(uR, uI, phiR, phiI, dynTimestep) ? 1 : 0;
					if (acceptNewParams != 1)
					{
						Log("PARAMETERS NOT ACCEPTED", ERROR);
					}
				}
				MPIMethods::BroadcastValue(&acceptNewParams);
				if (acceptNewParams != 1)
				{
					if (processRank == rootRank)
					{
						AllLocalEnergyR.pop_back();
						AllLocalEnergyI.pop_back();
						AllLocalOperators.pop_back();
						AllOtherExpectationValues.pop_back();
						AllParametersR.pop_back();
						AllParametersI.pop_back();
						uR = uRList.back();
						uI = uIList.back();
						phiR = phiRList.back();
						phiI = phiIList.back();
					}
					BroadcastNewParameters(uR, uI, &phiR, &phiI);
				}
			}
			else
			{
				acceptNewParams = 1;
			}
		}

		///////////////////////////////
		/// parameters are accepted ///
		/// write data to file      ///
		///////////////////////////////
		if (processRank == rootRank)
		{
			// Write current parameters
			if (WRITE_SINGLE_FILES == 1 && step % WRITE_EVERY_NTH_STEP_TO_FILE == 0)
			{
				WriteDataToFile(uR, "parametersR" + to_string(step), "parameterR, phiR=" + to_string(phiR) + ", wf=" + to_string(sys->GetWf()));
				WriteDataToFile(uI, "parametersI" + to_string(step), "parameterI, phiI=" + to_string(phiI));
			}
			// Write every nth expectation value calculated until now
			if (step % WRITE_EVERY_NTH_STEP_TO_FILE == 0)
			{
				//WriteDataToFile(AllLocalEnergyR, "AllLocalEnergyR", "ER", WRITE_EVERY_NTH_STEP_TO_FILE);
				//WriteDataToFile(AllLocalEnergyI, "AllLocalEnergyI", "EI", WRITE_EVERY_NTH_STEP_TO_FILE);
				//WriteDataToFile(AllLocalOperators, "AllLocalOperators", "<O_k>", WRITE_EVERY_NTH_STEP_TO_FILE);
				//WriteDataToFile(AllOtherExpectationValues, "AllOtherExpectationValues", "kinetic, potential, wf, g(r)", WRITE_EVERY_NTH_STEP_TO_FILE);
				//WriteDataToFile(AllParametersR, "AllParametersR", "uR", WRITE_EVERY_NTH_STEP_TO_FILE);
				//WriteDataToFile(AllParametersI, "AllParametersI", "uI", WRITE_EVERY_NTH_STEP_TO_FILE);

				AppendDataToFile(AllLocalEnergyR.back(), "AllLocalEnergyR");
				AppendDataToFile(AllLocalEnergyI.back(), "AllLocalEnergyI");
				AppendDataToFile(AllLocalOperators.back(), "AllLocalOperators");
				AppendDataToFile(AllOtherExpectationValues.back(), "AllOtherExpectationValues");
				AppendDataToFile(AllParametersR.back(), "AllParametersR");
				AppendDataToFile(AllParametersI.back(), "AllParametersI");
			}
		}

		// check if simulation should be cancelled or set back to a previous timestep
		int cancel = 0;
		int configChanged = 0;
		int setBack = 0;
		if (processRank == rootRank)
		{
			//cout << "#" << localEnergyR << " - " << std::fpclassify(localEnergyR) << endl;
			if (step % WRITE_EVERY_NTH_STEP_TO_FILE == 0 && FileExist("./stop")) //TODO: do not do this in every timestep. maybe every 30 seconds?
			{
				//TODO: read content of stop-file and interpret as new TOTALTIME instead of stopping immediately
				Log("Detected stop-file!", WARNING);
				cancel = 1;
			}
			else if (!std::isfinite(localEnergyR) || std::isnan(localEnergyR))
			{
				Log("Energy not finite", ERROR);
				cancel = 2;
			}
			else if (nrOfAcceptParameterTrials == maxNrOfAcceptParameterTrials)
			{
				Log("Parameters not accepted", ERROR);
				//cancel = 3;
				setBack = 5;
			}
			else if (step % WRITE_EVERY_NTH_STEP_TO_FILE == 0 && FileExist("./param")) //TODO: do not do this in every timestep. maybe every 30 seconds?
			{
				Log("Detected parameter change!", WARNING);
				PerformSystemParameterChange();
				configChanged = 1;
				RemoveFile("./param");
			}

			if (cancel != 0)
			{
				Log("Finishing simulation at t=" + to_string(currentTime));
			}
		}

		MPIMethods::BroadcastValue(&cancel);
		MPIMethods::BroadcastValue(&configChanged);
		MPIMethods::BroadcastValue(&setBack);
		if (cancel != 0)
		{
			TOTALTIME = currentTime;
			break; //break for (currentTime = 0; currentTime <= TOTALTIME; currentTime += TIMESTEP)
		}
		if (configChanged != 0)
		{
			BroadcastChangeableConfigItems();
			percentForExtraSampleToUpdate = modf(MC_NSTEPS * UPDATE_SAMPLES_PERCENT / 100.0, &nrOfSamplesToUpdate);
		}
		if (setBack != 0)
		{
			if (processRank == rootRank)
			{
				Log("SetBackNSteps", ERROR);
			}
			SetBackNSteps(setBack, R, uR, uI, phiR, phiI, step, dynTimestep);
			sys->CalculateWavefunction(R, uR, uI, phiR, phiI);
			for (int i = 0; i < MC_VERY_FIRST_NINITIALIZATIONSTEPS; i++)
			{
				DoMetropolisStep(R, uR, uI, phiR, phiI);
			}
			sys->CalculateWavefunction(R, uR, uI, phiR, phiI);
			setBack = 0;
		}

		if (processRank == rootRank)
		{
			t.stop();
			Log("duration for full timestep: " + to_string(t.duration()) + " ms");
		}
	}

	//if (processRank == rootRank)
	//{
	//	PrintData(uR);
	//	PrintData(uI);
	//}
	//cout << "phiR=" << phiR << ", phiI=" << phiI << endl;

	// Write data of whole simulation to files
	if (processRank == rootRank)
	{
		Log("Write last files ...");
		WriteDataToFile(AllLocalEnergyR, "AllLocalEnergyR", "ER");
		WriteDataToFile(AllLocalEnergyI, "AllLocalEnergyI", "EI");
		WriteDataToFile(AllLocalOperators, "AllLocalOperators", "<O_k>");
		WriteDataToFile(AllOtherExpectationValues, "AllOtherExpectationValues", "kinetic, potential, wf, g(r)");
		WriteDataToFile(AllParametersR, "AllParametersR", "uR");
		WriteDataToFile(AllParametersI, "AllParametersI", "uI");

		WriteDataToFile(AllLocalEnergyR, "AllLocalEnergyR_every100th", "ER", 100);
		WriteDataToFile(AllLocalEnergyI, "AllLocalEnergyI_every100th", "EI", 100);
		WriteDataToFile(AllLocalOperators, "AllLocalOperators_every100th", "ER", 100);
		WriteDataToFile(AllOtherExpectationValues, "AllOtherExpectationValues_every100th", "kinetic, potential, wf, g(r)", 100);
		WriteDataToFile(AllParametersR, "AllParametersR_every100th", "uR", 100);
		WriteDataToFile(AllParametersI, "AllParametersI_every100th", "uI", 100);

		WriteDataToFile(times, "AA_times", "t");
	}

	// Write random number gnerator states to file
	if (processRank == rootRank)
	{
		int tmp = 0;
		tmp = system(("mkdir " + OUT_DIR + "random").c_str());
		cout << tmp << endl;
	}
	MPIMethods::Barrier(); // wait for main process to create directory
	WriteRandomGeneratorStatesToFile("random/state");
	if (processRank == rootRank)
	{
		cout << "uniform: " << random01() << endl << "normal: " << randomNormal() << endl;
		cout << "uniform: " << random01() << endl << "normal: " << randomNormal() << endl;
		cout << "uniform: " << random01() << endl << "normal: " << randomNormal() << endl;
	}

	////////////////////////////////////////////////////
	/// additional calculations at end of simulation ///
	////////////////////////////////////////////////////

	nAcceptances = 0;
	nTrials = 0;
	if (USE_MEAN_FOR_FINAL_PARAMETERS == 1)
	{
		if (processRank == rootRank)
		{
			vector<double> finaluR = Mean(uRList);
			vector<double> finaluI = Mean(uIList);
			for (int i = 0; i < N_PARAM; i++)
			{
				uR[i] = finaluR[i];
				uI[i] = finaluI[i];
			}
			phiR = Mean(phiRList);
			phiI = Mean(phiIList);
			//WriteDataToFile(uRList, "finalParameterList", "uR[n]");
		}
	}
	BroadcastNewParameters(uR, uI, &phiR, &phiI);
	AlignCoordinates(R);
	ParallelCalculateAdditionalSystemProperties(R, uR, uI, phiR, phiI);
	if (processRank == rootRank)
	{
		WriteDataToFile(additionalSystemProperties, "AdditionalSystemProperties", "g(r), ...");
		WriteDataToFile(AllAdditionalSystemProperties, "AllAdditionalSystemProperties", "g(r), ...");
	}

	// Write config file for successive simulations
	AlignCoordinates(R);
	if (processRank == rootRank)
	{
		int tmp = 0;
		tmp = system(("mkdir " + OUT_DIR + "coords").c_str());
		cout << tmp << endl;
	}
	MPIMethods::Barrier();
	WriteParticleInputFile("coords/particleconfiguration_" + to_string(N) + "_" + to_string(DIM) + "D_" + to_string(processRank), R);
	if (processRank == rootRank)
	{
		MC_NSTEPS = mc_nsteps_original;
		WriteParticleInputFile("AAFinish_particleconfiguration_" + to_string(N) + "_" + to_string(DIM) + "D", R);
		cout << "phiR=" << phiR << endl;
		//cout << "log(additionalSystemProperties[2])=" << log(additionalSystemProperties[2]) << endl;
		//phiR = phiR - log(additionalSystemProperties[2]);
		cout << "phiR=" << phiR << endl;
		for (int i = 0; i < N_PARAM; i++)
		{
			PARAMS_REAL[i] = uR[i];
			PARAMS_IMAGINARY[i] = uI[i];
		}
		PARAM_PHIR = phiR;
		PARAM_PHII = phiI;
		WriteConfig("AAFinish_tVMC.config", uR, uI, phiR, phiI);

		if (FileExist("./stop"))
		{
			cout << "Delete stop-file!" << endl;
			RemoveFile("./stop");
		}
	}

	//Log("free memory ...");
	//TODO: how to properly delete the IPhysicalSystem pointer?
	//delete sys;
	for (int i = 0; i < MC_NSTEPS; i++)
	{
		delete correlatedSamplingData[i];
	}

	//Log("finalize ...");
	MPIMethods::Barrier();
	MPI_Finalize();

	//Log("done");

	return 0;
}

int startVMCSamplerMPI(int argc, char** argv)
{
	//////////////////////////
	/// Init MPI processes ///
	//////////////////////////

	char processName[80];
	int processNameLength;

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &processRank);
	MPI_Comm_size(MPI_COMM_WORLD, &numOfProcesses);
	MPI_Get_processor_name(processName, &processNameLength);

	MPIMethods::numOfProcesses = numOfProcesses;
	MPIMethods::processRank = processRank;
	MPIMethods::rootRank = rootRank;

	MPIMethods::GetCPUAllocation(false);

	////////////////////////////////////////
	/// Read and Broadcast configuration ///
	////////////////////////////////////////

	RegisterAllConfigItems();
	if (processRank == rootRank)
	{
		Log("Master process started...");

		if (FileExist(configFilePath))
		{
			ReadConfig(configFilePath);
			if (requiredConfigVersion != CONFIG_VERSION)
			{
				Log("Config file is out of date.\nrequired: \"" + requiredConfigVersion + "\"\nfound: \"" + CONFIG_VERSION + "\"", ERROR);
				return 1;
			}
			CreateOutputDirectory();
			//PrintConfig();
		}
		else
		{
			cout << "config file not found. path=\"" << configFilePath << "\"" << endl;
			MPI_Abort(MPI_COMM_WORLD, 0);
			return 1;
		}
	}
	configDirectory = configFilePath.substr(0, configFilePath.find_last_of("\\/")) + "/"; //TODO: restrict all file access to main process
	MPIMethods::BroadcastValue(&OUT_DIR, 300);
	//cout << "configDirectory=\"" << configDirectory << "\"" << endl;
	BroadcastConfig();
	//PrintConfig();

	//////////////////////////////////
	/// Initialize Physical System ///
	//////////////////////////////////

	if (!InitializePhysicalSystem())
	{
		return 1;
	}

	/////////////////////////////
	/// other Initializations ///
	/////////////////////////////

	Init();

	vector<vector<double> > R(N);
	vector<double> uR(N_PARAM);
	vector<double> uI(N_PARAM);
	double phiR;
	double phiI;
	PARAMS_REAL.resize(N_PARAM, 0.0);
	PARAMS_IMAGINARY.resize(N_PARAM, 0.0);

	for (unsigned int i = 0; i < R.size(); i++)
	{
		R[i].resize(DIM);
	}

	InitCoordinateConfiguration(R);
	if (processRank == rootRank)
	{
		WriteParticleInputFile("particleconfiguration", R);
	}
	if (processRank == rootRank)
	{
		for (int i = 0; i < N_PARAM; i++)
		{
			uR[i] = PARAMS_REAL[i];
			uI[i] = PARAMS_IMAGINARY[i];
		}
		phiR = PARAM_PHIR;
		phiI = PARAM_PHII;
		//PrintData(uR);
		//PrintData(uI);

		if (processRank == rootRank)		// && USE_ADJUST_PARAMETERS == true)
		{
			//AdjustParameters(uR, uI, &phiR, &phiI);
		}

		WriteDataToFile(uR, "parametersR0", "parameterR");
		WriteDataToFile(uI, "parametersI0", "parameterI");
	}

	sys->InitSystem();
	PostSystemInit();

	////////////////////////
	/// start simulation ///
	////////////////////////

	Timer t;
	if (processRank == rootRank)
	{
		t.start();
	}

	vector<vector<double> > energies;

	vector<double> startValues(uR);
	vector<double> steps = { 0.001, 0.01 };
	for (int i = 0; i < 60; i++)
	{
		for (int j = 0; j < 300; j++)
		{
			nAcceptances = 0;
			nTrials = 0;

			BroadcastNewParameters(uR, uI, &phiR, &phiI);

			AlignCoordinates(R);

			ParallelUpdateExpectationValues(R, uR, uI, phiR, phiI);

			int step = i;
			if (processRank == rootRank)
			{
				WriteDataToFile(otherExpectationValues, "otherExpectationValues" + to_string(step), "Ekin, Ekin_cor, Epot, Epot_corr, wf, g(r)_1, ..., g(r)_100");
			}

			if (processRank == rootRank)
			{
				energies.push_back(uR);
				energies[energies.size() - 1].push_back(localEnergyR / (double) N);
			}

			if (processRank == rootRank)
			{
				cout << i << ":" << JoinVector(uR) << endl;
				cout << "localEnergyR/N=" << localEnergyR / (double) N << " (" << otherExpectationValues[0] / (double) N << " + " << otherExpectationValues[1] / (double) N << " + " << endl;
				cout << t.interval() << "ms " << endl;
			}

			if (processRank == rootRank)
			{
				AppendDataToFile(energies.back(), "energies");
			}
			uR[0] += steps[0];
		}
		uR[0] = startValues[0];
		uR[1] += steps[1];
	}

	if (processRank == rootRank)
	{
		WriteDataToFile(energies, "energies", "a, energy");
	}

	MPIMethods::Barrier(); // wait for main process to create directory
	WriteRandomGeneratorStatesToFile("random/state");

	nAcceptances = 0;
	nTrials = 0;

	AlignCoordinates(R);
	ParallelCalculateAdditionalSystemProperties(R, uR, uI, phiR, phiI);
	if (processRank == rootRank)
	{
		WriteDataToFile(additionalSystemProperties, "AdditionalSystemProperties", "g(r), ...");
		WriteDataToFile(AllAdditionalSystemProperties, "AllAdditionalSystemProperties", "g(r), ...");
	}

	//Log("free memory ...");
	//TODO: how to properly delete the IPhysicalSystem pointer?
	//delete sys;
	for (int i = 0; i < MC_NSTEPS; i++)
	{
		delete correlatedSamplingData[i];
	}

	//Log("finalize ...");
	MPIMethods::Barrier();
	MPI_Finalize();

	//Log("done");

	return 0;
}

void startVMCSampler()
{
	int processRank = 0;
	VMCSampler::generator = mt19937_64(processRank + 1);

	LBOX = 4;
	LBOX_R = 1.0 / LBOX;
	N = 64;
	DIM = 3;
	N_PARAM = 100;
	MC_NSTEPS = 500;
	MC_NTHERMSTEPS = 50;
	MC_NINITIALIZATIONSTEPS = 10;
	int numberOfSplines = N_PARAM + 2;
	double halfLength = LBOX / 2.0;
	VMCSampler::nodePointSpacing = halfLength / (double) (numberOfSplines - 3.0);

	vector<vector<double> > energies;
	//double e1 = VMCSampler::GetEnergy();

	for (int i = 0; i < MC_NINITIALIZATIONSTEPS; i++)
	{
		VMCSampler::DoMetropolisStep();
	}
	for (int i = 0; i < MC_NSTEPS; i++)
	{
		for (int j = 0; j < MC_NTHERMSTEPS; j++)
		{
			VMCSampler::DoMetropolisStep();
		}
		double e = VMCSampler::GetEnergy();
		energies.push_back( { e, -VMCSampler::energy1 / (double) N, -VMCSampler::energy2 / (double) N, -VMCSampler::energy2_1 / (double) N, -VMCSampler::energy2_2 / (double) N });
		cout << JoinVector(energies.back()) << endl;
	}

	vector<double> energy = Mean(energies);

	cout << "E=" << JoinVector(energy) << endl;
}

int main(int argc, char **argv)
{
	int val = -1;
	//cout << to_string(argc) << "|" << argv[1] << "|" << to_string(strcmp(argv[1], "-mpitest")) << "|" << endl;
	//cout << "#################################################" << endl;
	//cout << "#################################################" << endl;
	//Log("start");

	//startVMCSampler();
	//return 0;

	Potentials::SaveAllPotentialValues("/itpstore/gartner/Output/TDVMC/asterix/");

	if (argc == 0) //INFO: started with pbs on mach
	{
		//cout << "Starting on MACH (argc=0)" << endl << endl;
		configFilePath = "/home/k3501/k354522/tVMC/bin/tVMC.config";
		val = mainMPI(argc, argv);
	}
	else if (argc == 1) //INFO: started without specifying config-file. used for local execution
	{
		//cout << "Starting (argc=1)" << endl << endl;
		//SYSTEM_TYPE = "BulkQT";
		//SYSTEM_TYPE = "BulkSplines";
		//SYSTEM_TYPE = "BulkSplinesScaled";
		//SYSTEM_TYPE = "HardSphereBosons";
		//SYSTEM_TYPE = "HardSphereBosonsExp";
		//SYSTEM_TYPE = "HeDrop";
		//SYSTEM_TYPE = "BosonCluster";
		//SYSTEM_TYPE = "BosonClusterWithLog";
		SYSTEM_TYPE = "BosonClusterWithLogParam";
		SYSTEM_TYPE = "BosonMixtureCluster";
		//SYSTEM_TYPE = "NUBosonsBulk";
		//SYSTEM_TYPE = "PBosonsBulk";
		//SYSTEM_TYPE = "BosonsBulkDamped";
		if (SYSTEM_TYPE == "HeDrop")
		{
			//configFilePath = "/home/gartner/Sources/TDVMC/config/drop_20.config";
			configFilePath = "/home/gartner/Sources/TDVMC/config/X4Drop.config";
		}
		else if (SYSTEM_TYPE == "HeBulk")
		{
			configFilePath = "/home/gartner/Sources/TDVMC/config/bulk_64.config";
		}
		else if (SYSTEM_TYPE == "BulkSplines")
		{
			configFilePath = "/home/gartner/Sources/TDVMC/config/bulkSplines.config";
		}
		else if (SYSTEM_TYPE == "BulkSplinesScaled")
		{
			configFilePath = "/home/gartner/Sources/TDVMC/config/bulkSplinesScaled.config";
		}
		else if (SYSTEM_TYPE == "BulkQT")
		{
			//configFilePath = "/home/gartner/Sources/TDVMC/config/bulkQT.config";
			configFilePath = "/home/gartner/Sources/TDVMC/config/test64.config";
		}
		else if (SYSTEM_TYPE == "GaussianWavepacket")
		{
			configFilePath = "/home/gartner/Sources/TDVMC/config/wavepacket.config";
		}
		else if (SYSTEM_TYPE == "BosonsBulk")
		{
			configFilePath = "/home/gartner/Sources/TDVMC/config/BosonsBulk1D.config";
		}
		else if (SYSTEM_TYPE == "BosonsBulkDamped")
		{
			configFilePath = "/home/gartner/Sources/TDVMC/config/BosonsBulkDamped.config";
		}
		else if (SYSTEM_TYPE == "BosonCluster")
		{
			configFilePath = "/home/gartner/Sources/TDVMC/config/BosonCluster.config";
		}
		else if (SYSTEM_TYPE == "BosonClusterWithLog")
		{
			configFilePath = "/home/gartner/Sources/TDVMC/config/BosonClusterWithLog.config";
		}
		else if (SYSTEM_TYPE == "BosonClusterWithLogParam")
		{
			configFilePath = "/home/gartner/Sources/TDVMC/config/BosonClusterWithLogParam.config";
		}
		else if (SYSTEM_TYPE == "BosonMixtureCluster")
		{
			configFilePath = "/home/gartner/Sources/TDVMC/config/BosonMixtureCluster.config";
		}
		else if (SYSTEM_TYPE == "HardSphereBosons")
		{
			configFilePath = "/home/gartner/Sources/TDVMC/config/HardSphereBosons.config";
		}
		else if (SYSTEM_TYPE == "HardSphereBosonsExp")
		{
			configFilePath = "/home/gartner/Sources/TDVMC/config/HardSphereBosonsExp.config";
		}
		else if (SYSTEM_TYPE == "NUBosonsBulk")
		{
			configFilePath = "/home/gartner/Sources/TDVMC/config/NUBosonsBulk1D.config";
			//configFilePath = "/home/gartner/Sources/TDVMC/config/NUBosonsBulk2D.config";
		}
		else if (SYSTEM_TYPE == "PBosonsBulk")
		{
			configFilePath = "/home/gartner/Sources/TDVMC/config/PBosonsBulk2D.config";
		}
		//if (processRank == rootRank) //INFO: no processRank assigned so far. this is done in mainMPI or startVMCSamplerMPI
		//{
		//	Log("configFilePath: " + configFilePath);
		//}
		val = mainMPI(argc, argv);
		//val = startVMCSamplerMPI(argc, argv);
	}
	else if (argc > 1)
	{
		//cout << "Starting (argc>1)" << endl << endl;
		if (strcmp(argv[1], "-test") == 0)
		{
			Test::CalculatePi(argc, argv);
			Test::TestVectorDisplacements();
			val = 0;
		}
		else if (strcmp(argv[1], "-pitest") == 0)
		{
			Test::CalculatePi(argc, argv);
			val = 0;
		}
		else if (strcmp(argv[1], "-vectortest") == 0)
		{
			Test::TestVectorDisplacements();
			val = 0;
		}
		else if (strcmp(argv[1], "-mpitest") == 0)
		{
			Test::TestMPI(argc, argv);
			val = 0;
		}
		else
		{
			//INFO: don't know why it does not work on zusie when configFilePath is passed as a third parameter to mainMPI -> therefore I use a global variable
			configFilePath = string(argv[1]);
			val = mainMPI(argc, argv);
			//val = startVMCSamplerMPI(argc, argv);
		}
	}

	//Log("exit");
	return val;
}

