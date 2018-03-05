#include <chrono>
#include "Constants.h"
#include <cstring>
#include <ctime>
#include <fstream>
#include "HeBulk.h"
#include "HeDrop.h"
#include <iomanip>
#include <json/json.h>
#include <math.h>
#include "MathOperators.h"

//#undef SEEK_SET
//#undef SEEK_END
//#undef SEEK_CUR

#include "mpi.h"
#include <signal.h>
#include <stdlib.h>
#include <string>
#include "test/PiCalculator.h"
#include "test/Tests.h"
#include <unistd.h>
#include "Utils.h"
#include <vector>

using namespace std;

IPhysicalSystem* sys;

bool USE_MEAN_FOR_FINAL_PARAMETERS = false;
//bool useNIC = false; //for drops
bool useNIC = true; //for bulk with PBC
//bool moveCOMToZero = false; //for drops with external potential
bool moveCOMToZero = true; //for drops without external potential

string OUTPUT_DIRECTORY;
string originalOutputDirectory;
string configDirectory;
int N;           	    		//number of particles
int DIM;     	        	 	//number of dimensions
int N_PARAM;	        	  	//number of parameters of trial function
double MC_STEP;
int MC_NSTEPS;
int MC_NTHERMSTEPS;
int MC_NINITIALIZATIONSTEPS;
int MC_NADDITIONALSTEPS;
int MC_NADDITIONALTHERMSTEPS;
int MC_NADDITIONALINITIALIZATIONSTEPS;
double RHO;
double RC;          			//cutoff for WF and LJ
double TIMESTEP;
double TOTALTIME;
int IMAGINARY_TIME;
double LBOX;
vector<double> PARAMS_REAL;
vector<double> PARAMS_IMAGINARY;
double PARAM_PHIR;
double PARAM_PHII;

int numOfProcesses = 1;
int rootRank = 0;
int processRank = 0;
int nAcceptances = 0;
int nTrials = 0;
double mc_nsteps;
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

vector<vector<double> > uRList;
vector<vector<double> > uIList;
vector<double> phiRList;
vector<double> phiIList;

vector<double> AllLocalEnergyR;
vector<vector<double> > AllOtherExpectationValues;
vector<vector<double> > AllParametersR;
vector<vector<double> > AllParametersI;
vector<vector<double> > AllAdditionalSystemProperties;

default_random_engine generator;
uniform_real_distribution<double> distUniform(0.0, 1.0);
normal_distribution<double> distNormal(0.0, 1.0);

void onSignalStop(int signum)
{
    if (signum == 10) //SIGUSR1
    {
        cout << "Received SIGUSR1!" << endl;
        cout << "Finishing simulation at t=" << currentTime << endl;
        TOTALTIME = currentTime;
    }
}

double random01()
{
	return distUniform(generator);
}

//double random01()
//{
//	return ((double) rand() / RAND_MAX);
//}

double randomNormal()
{
	return distNormal(generator);
}

double randomNormal(double sigma, double mu)
{
	return randomNormal() * sigma + mu;
}

void Log(string message)
{
	cout << "#" << setfill(' ') << setw(2) << processRank << "@" << get_cpu_id() << ": " << message << endl << flush;
}

void CreateOutputDirectory(string filePath)
{
	//INFO: append config options to output directory path, create the directory and copy the config file
	int tmp;
	originalOutputDirectory = OUTPUT_DIRECTORY;
	OUTPUT_DIRECTORY = OUTPUT_DIRECTORY + "step=" + to_string(MC_NSTEPS) + "_therm=" + to_string(MC_NTHERMSTEPS) + "_time=" + to_string(TIMESTEP) + "/";
	tmp = system(("rm -rf " + OUTPUT_DIRECTORY).c_str());
	tmp = system(("mkdir " + OUTPUT_DIRECTORY).c_str());
	CopyFile(filePath, OUTPUT_DIRECTORY + "vmc.config");
	cout << tmp << endl;
}

void PrintConfig()
{
	Log("===== Configuration ====");
	cout << "OUTPUT_DIRECTORY:" << OUTPUT_DIRECTORY << endl;
	cout << "N:" << N << endl;
	cout << "DIM:" << DIM << endl;
	cout << "N_PARAM:" << N_PARAM << endl;
	cout << "RHO:" << RHO << endl;
	cout << "RC:" << RC << endl;
	cout << "MC_STEP:" << MC_STEP << endl;
	cout << "MC_NSTEPS:" << MC_NSTEPS << endl;
	cout << "MC_NTHERMSTEPS:" << MC_NTHERMSTEPS << endl;
	cout << "MC_NINITIALIZATIONSTEPS:" << MC_NINITIALIZATIONSTEPS << endl;
	cout << "MC_NADDITIONALSTEPS:" << MC_NADDITIONALSTEPS << endl;
	cout << "MC_NADDITIONALTHERMSTEPS:" << MC_NADDITIONALTHERMSTEPS << endl;
	cout << "MC_NADDITIONALINITIALIZATIONSTEPS:" << MC_NADDITIONALINITIALIZATIONSTEPS << endl;
	cout << "TIMESTEP:" << TIMESTEP << endl;
	cout << "TOTALTIME:" << TOTALTIME << endl;
	cout << "IMAGINARY_TIME:" << IMAGINARY_TIME << endl;
	cout << "PARAMS_REAL" << JoinVector(PARAMS_REAL) << endl;
	cout << "PARAMS_IMAGINARY" << JoinVector(PARAMS_IMAGINARY) << endl;
	cout << "PARAM_PHIR:" << PARAM_PHIR << endl;
	cout << "PARAM_PHII:" << PARAM_PHII << endl;
	cout << "=============================" << endl << endl;
}

void ReadConfig(string filePath)
{
	cout << "file exits: " << FileExist(filePath) << "." << endl;
	Json::Value configData;
	Json::Reader configReader;
	ifstream configFile(filePath, ifstream::binary);
	configReader.parse(configFile, configData, false);

	OUTPUT_DIRECTORY = configData["OUTPUT_DIRECTORY"] == 0 ? OUTPUT_DIRECTORY : configData["OUTPUT_DIRECTORY"].asString();
	N = !configData["N"] ? N : configData["N"].asInt();
	DIM = !configData["DIM"] ? DIM : configData["DIM"].asInt();
	N_PARAM = !configData["N_PARAM"] ? N_PARAM : configData["N_PARAM"].asInt();
	RHO = !configData["RHO"] ? RHO : configData["RHO"].asDouble();
	RC = !configData["RC"] ? RC : configData["RC"].asDouble();
	MC_STEP = !configData["MC_STEP"] ? MC_STEP : configData["MC_STEP"].asDouble();
	MC_NSTEPS = !configData["MC_NSTEPS"] ? MC_NSTEPS : configData["MC_NSTEPS"].asInt();
	MC_NTHERMSTEPS = !configData["MC_NTHERMSTEPS"] ? MC_NTHERMSTEPS : configData["MC_NTHERMSTEPS"].asInt();
	MC_NINITIALIZATIONSTEPS = !configData["MC_NINITIALIZATIONSTEPS"] ? MC_NINITIALIZATIONSTEPS : configData["MC_NINITIALIZATIONSTEPS"].asInt();
	MC_NADDITIONALSTEPS = !configData["MC_NADDITIONALSTEPS"] ? MC_NADDITIONALSTEPS : configData["MC_NADDITIONALSTEPS"].asInt();
	MC_NADDITIONALTHERMSTEPS = !configData["MC_NADDITIONALTHERMSTEPS"] ? MC_NADDITIONALTHERMSTEPS : configData["MC_NADDITIONALTHERMSTEPS"].asInt();
	MC_NADDITIONALINITIALIZATIONSTEPS = !configData["MC_NADDITIONALINITIALIZATIONSTEPS"] ? MC_NADDITIONALINITIALIZATIONSTEPS : configData["MC_NADDITIONALINITIALIZATIONSTEPS"].asInt();
	TIMESTEP = !configData["TIMESTEP"] ? TIMESTEP : configData["TIMESTEP"].asDouble();
	TOTALTIME = !configData["TOTALTIME"] ? TOTALTIME : configData["TOTALTIME"].asDouble();
	IMAGINARY_TIME = !configData["IMAGINARY_TIME"] ? IMAGINARY_TIME : configData["IMAGINARY_TIME"].asInt();

	Json::Value paramsR = configData["PARAMS_REAL"];
	for (unsigned int i = 0; i < paramsR.size(); i++)
	{
		PARAMS_REAL.push_back(paramsR[i].asDouble());
	}
	Json::Value paramsI = configData["PARAMS_IMAGINARY"];
	for (unsigned int i = 0; i < paramsI.size(); i++)
	{
		PARAMS_IMAGINARY.push_back(paramsI[i].asDouble());
	}
	PARAM_PHIR = !configData["PARAM_PHIR"] ? PARAM_PHIR : configData["PARAM_PHIR"].asDouble();
	PARAM_PHII = !configData["PARAM_PHII"] ? PARAM_PHII : configData["PARAM_PHII"].asDouble();
}

void WriteConfig(string fileName, double* uR, double* uI, double phiR, double phiI)
{
	ofstream configFile;
	configFile.open(OUTPUT_DIRECTORY + fileName, ios::out);
	configFile.precision(8);

	configFile << "{" << endl;
	configFile << "\t" << "\"" << "OUTPUT_DIRECTORY" << "\"" << " : " << "\"" << originalOutputDirectory << "\"," << endl;
	configFile << "\t" << "\"" << "N" << "\"" << " : " << N << "," << endl;
	configFile << "\t" << "\"" << "DIM" << "\"" << " : " << DIM << "," << endl;
	configFile << "\t" << "\"" << "N_PARAM" << "\"" << " : " << N_PARAM << "," << endl;
	configFile << "\t" << "\"" << "RHO" << "\"" << " : " << RHO << "," << endl;
	configFile << "\t" << "\"" << "RC" << "\"" << " : " << RC << "," << endl;
	configFile << "\t" << "\"" << "MC_STEP" << "\"" << " : " << MC_STEP << "," << endl;
	configFile << "\t" << "\"" << "MC_NSTEPS" << "\"" << " : " << MC_NSTEPS << "," << endl;
	configFile << "\t" << "\"" << "MC_NTHERMSTEPS" << "\"" << " : " << MC_NTHERMSTEPS << "," << endl;
	configFile << "\t" << "\"" << "MC_NINITIALIZATIONSTEPS" << "\"" << " : " << MC_NINITIALIZATIONSTEPS << "," << endl;
	configFile << "\t" << "\"" << "MC_NADDITIONALSTEPS" << "\"" << " : " << MC_NADDITIONALSTEPS << "," << endl;
	configFile << "\t" << "\"" << "MC_NADDITIONALTHERMSTEPS" << "\"" << " : " << MC_NADDITIONALTHERMSTEPS << "," << endl;
	configFile << "\t" << "\"" << "MC_NADDITIONALINITIALIZATIONSTEPS" << "\"" << " : " << MC_NADDITIONALINITIALIZATIONSTEPS << "," << endl;
	configFile << "\t" << "\"" << "TIMESTEP" << "\"" << " : " << TIMESTEP << "," << endl;
	configFile << "\t" << "\"" << "TOTALTIME" << "\"" << " : " << TOTALTIME << "," << endl;
	configFile << "\t" << "\"" << "IMAGINARY_TIME" << "\"" << " : " << IMAGINARY_TIME << "," << endl;

	configFile << "\t" << "\"" << "PARAMS_REAL" << "\"" << " : " << endl;
	configFile << "\t" << "[" << endl;
	configFile << "\t\t";
	for (int i = 0; i < N_PARAM; i++)
	{
		configFile << fixed << uR[i];
		if (i != N_PARAM - 1)
		{
			configFile << ", ";
		}
	}
	configFile << endl;
	configFile << "\t" << "]," << endl;

	configFile << "\t" << "\"" << "PARAMS_IMAGINARY" << "\"" << " : " << endl;
	configFile << "\t" << "[" << endl;
	configFile << "\t\t";
	for (int i = 0; i < N_PARAM; i++)
	{
		configFile << fixed << uI[i];
		if (i != N_PARAM - 1)
		{
			configFile << ", ";
		}
	}
	configFile << endl;
	configFile << "\t" << "]," << endl;

	configFile << "\t" << "\"" << "PARAM_PHIR" << "\"" << " : " << phiR << "," << endl;
	configFile << "\t" << "\"" << "PARAM_PHII" << "\"" << " : " << phiI << endl;

	configFile << "}";

	configFile.close();
}

void WriteConfigJson(string fileName, double* uR, double* uI, double phiR, double phiI)
{
	ofstream configFile;
	Json::Value event;
    Json::Value paramsR(Json::arrayValue);
    Json::Value paramsI(Json::arrayValue);

    event["OUTPUT_DIRECTORY"]  = OUTPUT_DIRECTORY;
    event["N"]  = N;
    event["DIM"]  = DIM;
    event["N_PARAM"]  = N_PARAM;
    event["RHO"]  = RHO;
    event["RC"]  = RC;
    event["MC_STEP"]  = MC_STEP;
    event["MC_NSTEPS"]  = MC_NSTEPS;
    event["MC_NTHERMSTEPS"]  = MC_NTHERMSTEPS;
    event["MC_NINITIALIZATIONSTEPS"]  = MC_NINITIALIZATIONSTEPS;
    event["MC_NADDITIONALSTEPS"]  = MC_NADDITIONALSTEPS;
    event["MC_NADDITIONALTHERMSTEPS"]  = MC_NADDITIONALTHERMSTEPS;
    event["MC_NADDITIONALINITIALIZATIONSTEPS"]  = MC_NADDITIONALINITIALIZATIONSTEPS;
    event["TIMESTEP"]  = TIMESTEP;
    event["TOTALTIME"]  = TOTALTIME;
    event["IMAGINARY_TIME"]  = IMAGINARY_TIME;

    for (int i = 0; i < N_PARAM; i++)
    {
    	paramsR.append(Json::Value(uR[i]));
    	paramsI.append(Json::Value(uI[i]));
    }
    event["PARAMS_REAL"] = paramsR;
    event["PARAMS_IMAGINARY"] = paramsI;

    event["PARAM_PHIR"]  = phiR;
    event["PARAM_PHII"]  = phiI;

	configFile.open(OUTPUT_DIRECTORY + fileName, ios::out);
	configFile << event << endl;
	configFile.close();
}

void WriteParticleInputFile(string fileName, double** R)
{
	vector<vector<double> > particlePositionList = {{}};
	for (int i = 0; i < N; i++)
	{
		for (int a = 0; a < DIM; a++)
		{
			particlePositionList[0].push_back(R[i][a]);
		}
	}
	WriteDataToFile(particlePositionList, fileName, "", 1, false);
}

void ReadRandomGeneratorStatesFromFile(string fileNamePrefix)
{
	ifstream fGenerator(configDirectory + fileNamePrefix + "_generator.dat");
	fGenerator >> generator;
	fGenerator.close();
	ifstream fUniform(configDirectory + fileNamePrefix + "_uniform.dat");
	fUniform >> distUniform;
	fUniform.close();
	ifstream fNormal(configDirectory + fileNamePrefix + "_normal.dat");
	fNormal >> distNormal;
	fNormal.close();
}

void WriteRandomGeneratorStatesToFile(string fileNamePrefix)
{
	ofstream fGenerator(OUTPUT_DIRECTORY + fileNamePrefix + "_generator.dat");
	fGenerator << generator;
	fGenerator.close();
	ofstream fUniform(OUTPUT_DIRECTORY + fileNamePrefix + "_uniform.dat");
	fUniform << distUniform;
	fUniform.close();
	ofstream fNormal(OUTPUT_DIRECTORY + fileNamePrefix + "_normal.dat");
	fNormal << distNormal;
	fNormal.close();
}

void BroadcastConfig()
{
	char dir[200];
	strcpy(dir, OUTPUT_DIRECTORY.c_str());
	MPI_Bcast(dir, 200, MPI_CHAR, rootRank, MPI_COMM_WORLD);
	OUTPUT_DIRECTORY = string(dir);
	MPI_Bcast(&N, 1, MPI_INT, rootRank, MPI_COMM_WORLD);
	MPI_Bcast(&DIM, 1, MPI_INT, rootRank, MPI_COMM_WORLD);
	MPI_Bcast(&N_PARAM, 1, MPI_INT, rootRank, MPI_COMM_WORLD);
	MPI_Bcast(&RHO, 1, MPI_DOUBLE, rootRank, MPI_COMM_WORLD);
	MPI_Bcast(&RC, 1, MPI_DOUBLE, rootRank, MPI_COMM_WORLD);
	MPI_Bcast(&MC_STEP, 1, MPI_DOUBLE, rootRank, MPI_COMM_WORLD);
	MPI_Bcast(&MC_NSTEPS, 1, MPI_INT, rootRank, MPI_COMM_WORLD);
	MPI_Bcast(&MC_NTHERMSTEPS, 1, MPI_INT, rootRank, MPI_COMM_WORLD);
	MPI_Bcast(&MC_NINITIALIZATIONSTEPS, 1, MPI_INT, rootRank, MPI_COMM_WORLD);
	MPI_Bcast(&MC_NADDITIONALSTEPS, 1, MPI_INT, rootRank, MPI_COMM_WORLD);
	MPI_Bcast(&MC_NADDITIONALTHERMSTEPS, 1, MPI_INT, rootRank, MPI_COMM_WORLD);
	MPI_Bcast(&MC_NADDITIONALINITIALIZATIONSTEPS, 1, MPI_INT, rootRank, MPI_COMM_WORLD);
	MPI_Bcast(&TIMESTEP, 1, MPI_DOUBLE, rootRank, MPI_COMM_WORLD);
	MPI_Bcast(&TOTALTIME, 1, MPI_DOUBLE, rootRank, MPI_COMM_WORLD);
	MPI_Bcast(&IMAGINARY_TIME, 1, MPI_INT, rootRank, MPI_COMM_WORLD);
}

void BroadcastNewParameters(double* uR, double* uI, double* phiR, double* phiI)
{
	//TODO: check if phiR and phiI are broadcasted correctly
	MPI_Bcast(uR, N_PARAM, MPI_DOUBLE, rootRank, MPI_COMM_WORLD);
	MPI_Bcast(uI, N_PARAM, MPI_DOUBLE, rootRank, MPI_COMM_WORLD);
	MPI_Bcast(phiR, 1, MPI_DOUBLE, rootRank, MPI_COMM_WORLD);
	MPI_Bcast(phiI, 1, MPI_DOUBLE, rootRank, MPI_COMM_WORLD);
}

void BroadcastValue(double* value)
{
	MPI_Bcast(value, 1, MPI_DOUBLE, rootRank, MPI_COMM_WORLD);
}

void BroadcastValues(double* values, int count)
{
	MPI_Bcast(values, count, MPI_DOUBLE, rootRank, MPI_COMM_WORLD);
}

void BroadcastValues(vector<double>& values)
{
	double* valueArray = new double[values.size()];
	for (unsigned int i = 0; i < values.size(); i++)
	{
		valueArray[i] = values[i];
	}
	MPI_Bcast(valueArray, values.size(), MPI_DOUBLE, rootRank, MPI_COMM_WORLD);
	for (unsigned int i = 0; i < values.size(); i++)
	{
		values[i] = valueArray[i];
	}
	delete[] valueArray;
}

vector<double> ReduceToMinMaxMean(double* data)
{
	vector<double> values = {};
	double own = 0;
	double min = 0;
	double max = 0;
	double sum = 0;

	own = *data;
	MPI_Reduce(&own, &min, 1, MPI_DOUBLE, MPI_MIN, rootRank, MPI_COMM_WORLD);
	MPI_Reduce(&own, &max, 1, MPI_DOUBLE, MPI_MAX, rootRank, MPI_COMM_WORLD);
	MPI_Reduce(&own, &sum, 1, MPI_DOUBLE, MPI_SUM, rootRank, MPI_COMM_WORLD);
	if (processRank == rootRank)
	{
		values = { min, max, sum / (double)numOfProcesses };
	}
	return values;
}

void ReduceValues(double* data)
{
	double ownValues;
	double reducedValues;

	ownValues = *data;
	MPI_Reduce(&ownValues, &reducedValues, 1, MPI_DOUBLE, MPI_SUM, rootRank, MPI_COMM_WORLD);
	if (processRank == rootRank)
	{
		*data = reducedValues / (double)numOfProcesses;
	}
}

void ReduceValues(vector<double>& data)
{
	double* ownValues = new double[data.size()];
	double* reducedValues = new double[data.size()];

	for (unsigned int i = 0; i < data.size(); i++)
	{
		ownValues[i] = data[i];
	}
	MPI_Reduce(ownValues, reducedValues, data.size(), MPI_DOUBLE, MPI_SUM, rootRank, MPI_COMM_WORLD);
	//INFO: only root process holds the average values from all processes.
	if (processRank == rootRank)
	{
		for (unsigned int i = 0; i < data.size(); i++)
		{
			data[i] = reducedValues[i] / (double)numOfProcesses;
		}
	}
	delete[] ownValues;
	delete[] reducedValues;
}

void ReduceValues(vector<vector<double> >& data)
{
	for (unsigned int i = 0; i < data.size(); i++)
	{
		ReduceValues(data[i]);
	}
}

void Init()
{
	if (FileExist(configDirectory + "state" + "_generator.dat") && FileExist(configDirectory + "state" + "_uniform.dat") && FileExist(configDirectory + "state" + "_normal.dat"))
	{
		ReadRandomGeneratorStatesFromFile("state");
	}
	else
	{
		generator = default_random_engine(processRank + 1);
		srand(processRank + 1);
	}
	LBOX = pow((N / RHO), 1.0/DIM); //box with dimensions [-L/2, L/2]
	nAcceptances = 0;
	nTrials = 0;

    mc_nsteps = (double)MC_NSTEPS;
    mc_nadditionalsteps = (double)MC_NADDITIONALSTEPS;
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

void WriteParticlesToFile(double** R, string ending)
{
	WriteDataToFile(R, N, DIM, "position" + ending, "x, y, z");
}

bool LoadLastPositionsFromFile(string filename, double** R)
{
	bool successful = false;
	string line;
	string prevline;
	vector<string> coordinates;
	ifstream file;
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
		for (unsigned int i = 0; i < coordinates.size(); i += 3)
		{
			R[j][0] = stod(coordinates[i]);
			R[j][1] = stod(coordinates[i + 1]);
			R[j][2] = stod(coordinates[i + 2]);

			//double factor = 1.2;
			//R[j][0] *= factor;
			//R[j][1] *= factor;
			//R[j][2] *= factor;
			j++;
		}
		successful = true;
	}
	file.close();
	if (successful)
	{
		CopyFile("./" + filename + ".csv", OUTPUT_DIRECTORY + "particleconfiguration.csv");
	}
	return successful;
}

void InitCoordinateConfiguration(double** R)
{
	string filename = "particleconfiguration_" + to_string(N);
	if (!LoadLastPositionsFromFile(filename, R))
	{
		cout << "LBOX=" << LBOX << endl;
		int type = 2;
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
			double n = round(pow(N, 1.0/3.0));
			double l = pow(LBOX, 1.0/3.0);
			l *= 0.5;
			int p = 0;
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
		else if (type == 2)
		{
			//Drop
			Log("init coordinates for drop");
			int i = 0;
			double l = pow(LBOX, 1.0/3.0);
			vector<int> positions;
			while(i < N)
			{
				positions = {0, 0, 0};
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

void MoveCoordinatesToFirstCell(double** R)
{
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < DIM; j++)
		{
			R[i][j] = GetCoordinateNIC(R[i][j]);
		}
	}
}

void MoveCenterOfMassToZero(double** R)
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

void DoMetropolisStep(double** R, double* uR, double* uI, double phiR, double phiI)
{
	double p;
	bool sampleOkay = true;
	int randomParticle = ((int)floor(random01() * N)) % N;
	double* oldPosition = new double[DIM];
	for (int i = 0; i < DIM; i++)
	{
		oldPosition[i] = R[randomParticle][i];
		//R[randomParticle][i] += (random01() - 0.5) * MC_STEP;
		R[randomParticle][i] += randomNormal(MC_STEP, 0.0);
	}
	double wfQuotient = sys->GetWFQuotient(R, uR, uI, phiR, phiI, randomParticle, oldPosition);

	if (wfQuotient == 0 || sys->GetWfNew() == 0 || sys->GetWf() == 0)
	{
		//cout << "something == 0 !!!!!!!!!!!!!!!!!!!!!!!!" << endl;
		//cout << "wfQuotient=" << wfQuotient << ", wfNew=" << Sys::wfNew << ", wf=" << Sys::wf << ", nTrials=" << nTrials << endl;
		//sleep(1);
	}
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
	delete[] oldPosition;
}

void UpdateExpectationValues(double** R, double* uR, double* uI, double phiR, double phiI, bool intermediateStep = false)
{
	//vector<vector<double> > singlelocalOperators; // for all O_k
	vector<double> singlelocalEnergyR; // for all E^R
	//vector<double> singlelocalEnergyI; // for all E^I
	//vector<vector<vector<double> > > singlelocalOperatorsMatrix; // for all O_k O_j
	//vector<vector<double> > singlelocalOperatorlocalEnergyR; // for all O_k E^R
	//vector<vector<double> > singlelocalOperatorlocalEnergyI; // for all O_k E^I

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

		//INFO: consumes too much memory
		//singlelocalOperators.push_back(Sys::localOperators);
		singlelocalEnergyR.push_back(sys->GetLocalEnergyR());
		//singlelocalEnergyI.push_back(Sys::localEnergyI);
		//singlelocalOperatorsMatrix.push_back(Sys::localOperatorsMatrix);
		//singlelocalOperatorlocalEnergyR.push_back(Sys::localOperatorlocalEnergyR);
		//singlelocalOperatorlocalEnergyI.push_back(Sys::localOperatorlocalEnergyI);

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
			if (processRank == rootRank)
			{
				cout << "." << flush;
			}
		}
	}
	if (processRank == rootRank)
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
}

void ParallelUpdateExpectationValues(double** R, double* uR, double* uI, double phiR, double phiI, bool intermediateStep = false)
{
	chrono::high_resolution_clock::time_point t1;
	chrono::high_resolution_clock::time_point t2;
	int duration;
	double dblDuration;
	t1 = chrono::high_resolution_clock::now();

	UpdateExpectationValues(R, uR, uI, phiR, phiI, intermediateStep);

	t2 = chrono::high_resolution_clock::now();
	duration = chrono::duration_cast<chrono::milliseconds>( t2 - t1 ).count();
	dblDuration = (double) duration;
	vector<double> timings = ReduceToMinMaxMean(&dblDuration);
	if (processRank == rootRank && !intermediateStep)
	{
		Log("duration: min = " + to_string(timings[0]) + " ms");
		Log("          max = " + to_string(timings[1]) + " ms");
		Log("          <t> = " + to_string(timings[2]) + " ms");
	}

	ReduceValues(localOperators);
	ReduceValues(&localEnergyR);
	ReduceValues(&localEnergyI);
	ReduceValues(localOperatorsMatrix);
	ReduceValues(localOperatorlocalEnergyR);
	ReduceValues(localOperatorlocalEnergyI);
	ReduceValues(otherExpectationValues);
}

void PrintParameters(double* u)
{
	for (int i = 0; i < N_PARAM; i++)
	{
		cout << u[i] << ", ";
	}
	cout << endl;
}

void NormalizeParameters(double* uR, double* uI, double *phiR, double *phiI)
{
	//double sum = 0;
	//double factor = 1;
	//for (int i = 0; i < N_PARAM; i++)
	//{
	//	sum += u[i];
	//}
	//factor = N_PARAM / (2.0 * sum);
	//for (int i = 0; i < N_PARAM; i++)
	//{
	//	u[i] *= factor;
	//}
	double factor = 2.0 / (N * (N - 1.0));
	for (int i = 0; i < N_PARAM; i++)
	{
		uR[i] *= factor;
	}
}

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
			energiesImag[i] = - localOperatorlocalEnergyR[i];
		}
		else
		{
			energiesReal[i] = - localOperatorlocalEnergyR[i];
			energiesImag[i] = - localOperatorlocalEnergyI[i];
		}
		for (int j = 0; j <= i; j++)
		{
			matrix[i][j] = localOperatorsMatrix[i][j];
			matrix[j][i] = matrix[i][j];
		}
	}
}

void BuildSystemOfEquationsForParameters(vector<vector<double> >& matrix, vector<double>& energiesReal, vector<double>& energiesImag)
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
			energiesImag[i] = - localOperatorlocalEnergyR[i] + localEnergyR * localOperators[i];
		}
		else
		{
			energiesReal[i] = - localOperatorlocalEnergyR[i] + localEnergyR * localOperators[i];
			energiesImag[i] = - localOperatorlocalEnergyI[i];
		}
		for (int j = 0; j <= i; j++)
		{
			matrix[i][j] = localOperatorsMatrix[i][j] - localOperators[i] * localOperators[j];
			matrix[j][i] = matrix[i][j];
		}
	}
}

void PerformCholeskyDecomposition(vector<vector<double> >& matrix)
{
	double sum = 0;
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
				cout << "!!!!! NOT POSITIVE SEMI DEFINITE !!!! @ i=" << i << ", j=" << j << endl;
			}
		}
	}
}

void SolveForParametersDot(vector<vector<double> >& matrix, vector<double>& energiesReal, vector<double>& energiesImag, vector<double>& resultReal, vector<double>& resultImag)
{
	double sumReal = 0;
	double sumImag = 0;
	vector<double> tmpReal;
	vector<double> tmpImag;
	tmpReal.resize(N_PARAM);
	tmpImag.resize(N_PARAM);

	resultReal.resize(N_PARAM);
	resultImag.resize(N_PARAM);

	for (int i = 0; i < N_PARAM; i++)
	{
		sumReal = 0;
		sumImag = 0;
		for (int j = 0; j < i; j++)
		{
			sumReal += matrix[i][j] * tmpReal[j];
			sumImag += matrix[i][j] * tmpImag[j];
		}
		tmpReal[i] = 1.0 / matrix[i][i] * (energiesReal[i] - sumReal);
		tmpImag[i] = 1.0 / matrix[i][i] * (energiesImag[i] - sumImag);
	}

	for (int i = N_PARAM - 1; i >= 0; i--)
	{
		sumReal = 0;
		sumImag = 0;
		for (int j = N_PARAM - 1; j > i; j--)
		{
			sumReal += matrix[j][i] * resultReal[j];
			sumImag += matrix[j][i] * resultImag[j];
		}
		resultReal[i] = 1.0 / matrix[i][i] * (tmpReal[i] - sumReal);
		resultImag[i] = 1.0 / matrix[i][i] * (tmpImag[i] - sumImag);
		//cout << resultReal[i] << endl;
	}
}

void SolveForParametersDot(vector<vector<double> >& matrix, vector<double>& energiesReal, vector<double>& energiesImag, vector<double>& resultReal, vector<double>& resultImag, double *resultPhiReal, double *resultPhiImag)
{
	double sumReal = 0;
	double sumImag = 0;
	vector<double> tmpReal;
	vector<double> tmpImag;
	tmpReal.resize(N_PARAM);
	tmpImag.resize(N_PARAM);

	resultReal.resize(N_PARAM);
	resultImag.resize(N_PARAM);
	*resultPhiReal = 0;
	*resultPhiImag = 0;

	for (int i = 0; i < N_PARAM; i++)
	{
		sumReal = 0;
		sumImag = 0;
		for (int j = 0; j < i; j++)
		{
			sumReal += matrix[i][j] * tmpReal[j];
			sumImag += matrix[i][j] * tmpImag[j];
		}
		tmpReal[i] = 1.0 / matrix[i][i] * (energiesReal[i] - sumReal);
		tmpImag[i] = 1.0 / matrix[i][i] * (energiesImag[i] - sumImag);
	}

	for (int i = N_PARAM - 1; i >= 0; i--)
	{
		sumReal = 0;
		sumImag = 0;
		for (int j = N_PARAM - 1; j > i; j--)
		{
			sumReal += matrix[j][i] * resultReal[j];
			sumImag += matrix[j][i] * resultImag[j];
		}
		resultReal[i] = 1.0 / matrix[i][i] * (tmpReal[i] - sumReal);
		resultImag[i] = 1.0 / matrix[i][i] * (tmpImag[i] - sumImag);
		*resultPhiReal -= localOperators[i] * resultReal[i];
		*resultPhiImag -= localOperators[i] * resultImag[i];
		//cout << resultReal[i] << endl;
	}
	if (IMAGINARY_TIME == 0)
	{
		*resultPhiImag -= localEnergyR;
	}
	else
	{
		*resultPhiReal -= localEnergyR;
	}
}

void CalculateNextParametersEuler(double* uR, double* uI)
{
	if (processRank == rootRank)
	{
		vector<double> energiesReal; //rhs of the equation that corresponds to uDotR;
		vector<double> energiesImag; //rhs of the equation that corresponds to uDotI;
		vector<vector<double> > matrix;
		vector<double> uDotR;
		vector<double> uDotI;

		BuildSystemOfEquationsForParameters(matrix, energiesReal, energiesImag);

		//WriteDataToFile(matrix, "BBmatrix", "BBmatrix");
		//WriteDataToFile(energiesReal, "BBenergiesReal", "BBenergiesReal");
		//WriteDataToFile(energiesImag, "BBenergiesImag", "BBenergiesImag");

		PerformCholeskyDecomposition(matrix);
		SolveForParametersDot(matrix, energiesReal, energiesImag, uDotR, uDotI);

		//WriteDataToFile(matrix, "BBcholesky", "BBcholesky");
		//WriteDataToFile(uDotR, "BBuDotR", "BBuDotR");

		for (int i = 0; i < N_PARAM; i++)
		{
			uR[i] = uR[i] + uDotR[i] * TIMESTEP;
			//uI[i] = uI[i] + uDotI[i] * TIMESTEP;
		}
	}
}

void CalculateNextParametersEuler(double* uR, double* uI, double *phiR, double *phiI)
{
	if (processRank == rootRank)
	{
		vector<double> energiesReal; //rhs of the equation that corresponds to uDotR;
		vector<double> energiesImag; //rhs of the equation that corresponds to uDotI;
		vector<vector<double> > matrix;
		vector<double> uDotR;
		vector<double> uDotI;
		double phiDotR = 0;
		double phiDotI = 0;

		BuildSystemOfEquationsForParameters(matrix, energiesReal, energiesImag);

		//WriteDataToFile(matrix, "BBmatrix", "BBmatrix");
		//WriteDataToFile(energiesReal, "BBenergiesReal", "BBenergiesReal");
		//WriteDataToFile(energiesImag, "BBenergiesImag", "BBenergiesImag");

		PerformCholeskyDecomposition(matrix);
		SolveForParametersDot(matrix, energiesReal, energiesImag, uDotR, uDotI, &phiDotR, &phiDotI);

		//WriteDataToFile(matrix, "BBcholesky", "BBcholesky");
		//WriteDataToFile(uDotR, "BBuDotR", "BBuDotR");

		for (int i = 0; i < N_PARAM; i++)
		{
			uR[i] = uR[i] + uDotR[i] * TIMESTEP;
			uI[i] = uI[i] + uDotI[i] * TIMESTEP;
		}
		*phiR = *phiR + phiDotR * TIMESTEP;
		*phiI = *phiI + phiDotI * TIMESTEP;
	}
}

void CalculateNextParametersPC(double** R, double* uR, double* uI, double *phiR, double *phiI)
{
	vector<double> energiesReal; //rhs of the equation that corresponds to uDotR;
	vector<double> energiesImag; //rhs of the equation that corresponds to uDotI;
	vector<vector<double> > matrix;
	vector<double> uDotR;
	vector<double> uDotI;
	double phiDotR = 0;
	double phiDotI = 0;

	int PCsteps = 1;
	double* tmpUR;
	double* tmpUI;
	double tmpPhiR = 0;
	double tmpPhiI = 0;
	vector<double> nextUDotR;
	vector<double> nextUDotI;
	double nextPhiDotR = 0;
	double nextPhiDotI = 0;

	tmpUR = new double[N_PARAM];
	tmpUI = new double[N_PARAM];
	for (int i = 0; i < N_PARAM; i++)
	{
		tmpUR[i] = 0;
		tmpUI[i] = 0;
	}

	if (processRank == rootRank)
	{
		BuildSystemOfEquationsForParameters(matrix, energiesReal, energiesImag);
		PerformCholeskyDecomposition(matrix);
		SolveForParametersDot(matrix, energiesReal, energiesImag, uDotR, uDotI, &phiDotR, &phiDotI);
		for (int i = 0; i < N_PARAM; i++)
		{
			tmpUR[i] = uR[i] + uDotR[i] * TIMESTEP;
			//tmpUI[i] = uI[i] + uDotI[i] * TIMESTEP;
		}
		tmpPhiR = *phiR + phiDotR * TIMESTEP;
		//tmpPhiI = *phiI + phiDotI * TIMESTEP;
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
			BuildSystemOfEquationsForParameters(matrix, energiesReal, energiesImag);
			PerformCholeskyDecomposition(matrix);
			SolveForParametersDot(matrix, energiesReal, energiesImag, nextUDotR, nextUDotI, &nextPhiDotR, &nextPhiDotI);
			for (int i = 0; i < N_PARAM; i++)
			{
				tmpUR[i] = uR[i] + (uDotR[i] + nextUDotR[i]) * TIMESTEP / 2.0;
				//tmpUI[i] = uI[i] + (uDotI[i] + nextUDotI[i]) * TIMESTEP / 2.0;
			}
			tmpPhiR = *phiR + (phiDotR + nextPhiDotR) * TIMESTEP / 2.0;
			//tmpPhiI = *phiI + (phiDotI + nextPhiDotI) * TIMESTEP / 2.0;
		}
	}

	for (int i = 0; i < N_PARAM; i++)
	{
		uR[i] = tmpUR[i];
		//uI[i] = tmpUI[i];
	}
	*phiR = tmpPhiR;
	//*phiI = tmpPhiI;

	delete[] tmpUR;
	delete[] tmpUI;
}

void CalculateNextParametersRK4(double** R, double* uR, double* uI, double *phiR, double *phiI)
{
	vector<double> energiesReal; //rhs of the equation that corresponds to uDotR;
	vector<double> energiesImag; //rhs of the equation that corresponds to uDotI;
	vector<vector<double> > matrix;

	double* tmpUR;
	double* tmpUI;
	double tmpPhiR = 0;
	double tmpPhiI = 0;

	vector<vector<double> > uDotR;
	vector<vector<double> > uDotI;
	vector<double> phiDotR;
	vector<double> phiDotI;

	tmpUR = new double[N_PARAM];
	tmpUI = new double[N_PARAM];
    for (int i = 0; i < N_PARAM; i++)
    {
            tmpUR[i] = 0;
            tmpUI[i] = 0;
    }
	uDotR.resize(4);
	uDotI.resize(4);
	phiDotR.resize(4);
	phiDotI.resize(4);

	if (processRank == rootRank)
	{
		BuildSystemOfEquationsForParameters(matrix, energiesReal, energiesImag);
		PerformCholeskyDecomposition(matrix);
		SolveForParametersDot(matrix, energiesReal, energiesImag, uDotR[0], uDotI[0], &(phiDotR[0]), &(phiDotI[0]));
		for (int i = 0; i < N_PARAM; i++)
		{
			tmpUR[i] = uR[i] + uDotR[0][i] * TIMESTEP / 2.0;
			//tmpUI[i] = uI[i] + uDotI[0][i] * TIMESTEP / 2.0;
		}
		tmpPhiR = *phiR + phiDotR[0] * TIMESTEP / 2.0;
		//tmpPhiI = *phiI + phiDotI[0] * TIMESTEP / 2.0;
	}
	BroadcastNewParameters(tmpUR, tmpUI, &tmpPhiR, &tmpPhiI);
	ParallelUpdateExpectationValues(R, tmpUR, tmpUI, tmpPhiR, tmpPhiI, true);

	if (processRank == rootRank)
	{
		BuildSystemOfEquationsForParameters(matrix, energiesReal, energiesImag);
		PerformCholeskyDecomposition(matrix);
		SolveForParametersDot(matrix, energiesReal, energiesImag, uDotR[1], uDotI[1], &(phiDotR[1]), &(phiDotI[1]));
		for (int i = 0; i < N_PARAM; i++)
		{
			tmpUR[i] = uR[i] + uDotR[1][i] * TIMESTEP / 2.0;
			//tmpUI[i] = uI[i] + uDotI[1][i] * TIMESTEP / 2.0;
		}
		tmpPhiR = *phiR + phiDotR[1] * TIMESTEP / 2.0;
		//tmpPhiI = *phiI + phiDotI[1] * TIMESTEP / 2.0;
	}
	BroadcastNewParameters(tmpUR, tmpUI, &tmpPhiR, &tmpPhiI);
	ParallelUpdateExpectationValues(R, tmpUR, tmpUI, tmpPhiR, tmpPhiI, true);

	if (processRank == rootRank)
	{
		BuildSystemOfEquationsForParameters(matrix, energiesReal, energiesImag);
		PerformCholeskyDecomposition(matrix);
		SolveForParametersDot(matrix, energiesReal, energiesImag, uDotR[2], uDotI[2], &(phiDotR[2]), &(phiDotI[2]));
		for (int i = 0; i < N_PARAM; i++)
		{
			tmpUR[i] = uR[i] + uDotR[2][i] * TIMESTEP;
			//tmpUI[i] = uI[i] + uDotI[2][i] * TIMESTEP;
		}
		tmpPhiR = *phiR + phiDotR[2] * TIMESTEP;
		//tmpPhiI = *phiI + phiDotI[2] * TIMESTEP;
	}
	BroadcastNewParameters(tmpUR, tmpUI, &tmpPhiR, &tmpPhiI);
	ParallelUpdateExpectationValues(R, tmpUR, tmpUI, tmpPhiR, tmpPhiI, true);

	if (processRank == rootRank)
	{
		BuildSystemOfEquationsForParameters(matrix, energiesReal, energiesImag);
		PerformCholeskyDecomposition(matrix);
		SolveForParametersDot(matrix, energiesReal, energiesImag, uDotR[3], uDotI[3], &(phiDotR[3]), &(phiDotI[3]));
		for (int i = 0; i < N_PARAM; i++)
		{
			uR[i] = uR[i] + ((uDotR[0][i] + 2.0 * uDotR[1][i] + 2.0 * uDotR[2][i] + uDotR[3][i]) / 6.0) * TIMESTEP;
			//uI[i] = uI[i] + ((uDotI[0][i] + 2.0 * uDotI[1][i] + 2.0 * uDotI[2][i] + uDotI[3][i]) / 6.0) * TIMESTEP;
		}
		*phiR = *phiR + ((phiDotR[0] + 2.0 * phiDotR[1] + 2.0 * phiDotR[2] + phiDotR[3]) / 6.0) * TIMESTEP;
		//phiI = *phiI + ((phiDotI[0] + 2.0 * phiDotI[1] + 2.0 * phiDotI[2] + phiDotI[3]) / 6.0) * TIMESTEP;
	}

	delete[] tmpUR;
	delete[] tmpUI;
}

void CalculateNextParameters(double** R, double* uR, double* uI)
{
	int type = 0;
	if (type == 0)
	{
		CalculateNextParametersEuler(uR, uI);
	}
}

void CalculateNextParameters(double** R, double* uR, double* uI, double *phiR, double *phiI)
{
	chrono::high_resolution_clock::time_point t1;
	chrono::high_resolution_clock::time_point t2;
	int duration;
	if (processRank == rootRank)
	{
		t1 = chrono::high_resolution_clock::now();
	}

	int type = 0;
	if (type == 0)
	{
		CalculateNextParametersEuler(uR, uI, phiR, phiI);
	}
	else if (type == 1)
	{
		double** Rcopy;
		Rcopy = new double*[N];
		for (int i = 0; i < N; i++)
		{
			Rcopy[i] = new double[DIM];
		}
		Copy2DArray(R, Rcopy, N, DIM);
		CalculateNextParametersPC(Rcopy, uR, uI, phiR, phiI);
		for (int i = 0; i < N; i++)
		{
			delete[] Rcopy[i];
		}
		delete[] Rcopy;
	}
	else if (type == 2)
	{
		double** Rcopy;
		Rcopy = new double*[N];
		for (int i = 0; i < N; i++)
		{
			Rcopy[i] = new double[DIM];
		}
		Copy2DArray(R, Rcopy, N, DIM);
		CalculateNextParametersRK4(Rcopy, uR, uI, phiR, phiI);
		for (int i = 0; i < N; i++)
		{
			delete[] Rcopy[i];
		}
		delete[] Rcopy;
	}

	if (processRank == rootRank)
	{
		t2 = chrono::high_resolution_clock::now();
		duration = chrono::duration_cast<chrono::milliseconds>( t2 - t1 ).count();
		Log("DGL duration = " + to_string(duration) + " ms");
	}
}

void NormalizeWavefunction(double wf, double *phiR)
{
	*phiR = *phiR - log(wf);
}

void CalculateAdditionalSystemProperties(double** R, double* uR, double* uI, double phiR, double phiI)
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

void ParallelCalculateAdditionalSystemProperties(double** R, double* uR, double* uI, double phiR, double phiI)
{
	CalculateAdditionalSystemProperties(R, uR, uI, phiR, phiI);

	ReduceValues(additionalSystemProperties);
}

void AlignCoordinates(double** R)
{
	if (useNIC)
	{
		MoveCoordinatesToFirstCell(R);
	}
	else
	{
		if (moveCOMToZero)
		{
			MoveCenterOfMassToZero(R);
		}
	}
}

int mainSerial(int argc, char** argv)
{
	if (argc < 2)
	{
		cout << "no config file specified... argc=" << argc << endl;
		return 1;
	}
	ReadConfig(argv[1]);
	CreateOutputDirectory(argv[1]);
	//PrintConfig();
	Init();

	double** R;
	double** Rcopy;
	double* uR;
	double* uI;
	double phiR;
	double phiI;

	R = new double*[N];
	Rcopy = new double*[N];
	for (int i = 0; i < N; i++)
	{
		R[i] = new double[DIM];
		Rcopy[i] = new double[DIM];
	}
	uR = new double[N_PARAM];
	uI = new double[N_PARAM];

	InitCoordinateConfiguration(R);
	for (int i = 0; i < N_PARAM; i++)
	{
		//u[i] = PARAM_SPACE[i][0];
		//u[i] = 1 + sin(i / (N_PARAM / 2.0 / M_PI));
		//u[i] = 1.0;
		//u[i] = i / (double)N_PARAM - 0.2;
		//u[i] = random01();
		//u[i] = -1000.0 * exp(-(i+0.1)/4.0) + 50.0;
		//u[i] = exp(-100.0/pow(i, 2.0));
		//u[i] = pow(exp(-10.0/(double)i), 5);
		//u[i] = exp(-10000000.0/pow(i,5))*10.0;
		//u[i] = exp(-pow((60.0/(double)i), 5))*1000.0;
		//u[i] = exp(-pow(i/10.0, 5))*100;
		//if (i > N_PARAM / 2.0)
		//{
		//	u[i] = -10.0;
		//}
		//else
		//{
		//	u[i] = 1.0;
		//}
		//double ii = 2.7 * (i + N_PARAM / 3.0) / (double)N_PARAM;
		//u[i] = pow(ii, 2) * cos(ii) + 1.0;
		//double ii = i * M_PI / N_PARAM;
		//u[i] = sin(ii + 0.1) * exp(-5.0 * pow(ii, 3));
		//double ii = (i + 1.0) / (double)(N_PARAM / 1.3);
		//uR[i] = min(-1.0 / pow(ii, 1) + 1.0, 0.0);
		//uI[i] = 0.0;
		//uR[i] = 1.0;
		//uR[i] = 1.0;
		//uR[i] = 1.0 / (double)(i + 1.0);
		//uI[i] = 0.0;
		uR[i] = PARAMS_REAL[i];
		uI[i] = PARAMS_IMAGINARY[i];
	}
	phiR = PARAM_PHIR;
	phiI = PARAM_PHII;
	//phiR = -1.0;
	//phiI = 0.0;
	//NormalizeParameters(uR, uI, &phiR, &phiI);
	//PrintParameters(uR);
	//PrintParameters(uI);
	WriteDataToFile(uR, N_PARAM, "parametersR0", "parameterR");
	WriteDataToFile(uI, N_PARAM, "parametersI0", "parameterI");
	//exit(0);

	int step = 0;
	cout << "LBOX=" << LBOX << endl;
	sys->InitSystem();
	for (double t = 0; t <= TOTALTIME; t += TIMESTEP)
	{
		nAcceptances = 0;
		nTrials = 0;
		step++;

		UpdateExpectationValues(R, uR, uI, phiR, phiI);
		WriteDataToFile(localOperators, "AAlocalOperators", "localOperators");
		WriteDataToFile(localEnergyR, "AAlocalEnergyR", "localEnergyR");
		WriteDataToFile(localEnergyI, "AAlocalEnergyI", "localEnergyI");
		WriteDataToFile(localOperatorsMatrix, "AAlocalOperatorsMatrix", "localOperatorsMatrix");
		WriteDataToFile(localOperatorlocalEnergyR, "AAlocalOperatorlocalEnergyR", "localOperatorlocalEnergyR");
		WriteDataToFile(localOperatorlocalEnergyI, "AAlocalOperatorlocalEnergyI", "localOperatorlocalEnergyI");
		WriteDataToFile(otherExpectationValues, "AAotherExpectationValues", "otherExpectationValues");
		cout << "t=" << t << endl;
		cout << "localEnergyR=" << localEnergyR << endl;
		cout << "kineticEnergy=" << otherExpectationValues[0] << endl;
		cout << "potentialEnergy=" << otherExpectationValues[1] << endl;
		cout << "==========================================" << endl << endl;
		CalculateNextParameters(R, uR, uI, &phiR, &phiI);
		WriteDataToFile(uR, N_PARAM, "parametersR" + to_string(step), "parameterR");
		WriteDataToFile(uI, N_PARAM, "parametersI" + to_string(step), "parameterI");
	}

	for (int i = 0; i < N; i++)
	{
		delete[] R[i];
		delete[] Rcopy[i];
	}
	delete[] R;
	delete[] Rcopy;
	delete[] uR;
	delete[] uI;

	return 0;
}

int mainMPI(int argc, char** argv, string configFilePath)
{
	char processName[80];
	int processNameLength;

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &processRank);
	MPI_Comm_size(MPI_COMM_WORLD, &numOfProcesses);
	MPI_Get_processor_name(processName, &processNameLength);

	Log("running on cpu " + to_string(get_cpu_id()));

	if (processRank == rootRank)
	{
		Log("Master process started...");

		if (FileExist(configFilePath))
		{
			ReadConfig(configFilePath);
			CreateOutputDirectory(configFilePath);
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
	//cout << "configDirectory=\"" << configDirectory << "\"" << endl;
	BroadcastConfig();
	//PrintConfig();

	sys = new HeDrop(configDirectory);
	//sys = new HeBulk(configDirectory);
	Init();
	cout << "uniform: " << random01() << endl << "normal: " << randomNormal() << endl;
	cout << "uniform: " << random01() << endl << "normal: " << randomNormal() << endl;
	cout << "uniform: " << random01() << endl << "normal: " << randomNormal() << endl;

	double** R;
	double* uR;
	double* uI;
	double phiR;
	double phiI;

	R = new double*[N];
	for (int i = 0; i < N; i++)
	{
		R[i] = new double[DIM];
	}
	uR = new double[N_PARAM];
	uI = new double[N_PARAM];

	InitCoordinateConfiguration(R); //TODO: only init files in main process (at least when read from file)
	if (processRank == rootRank)
	{
		for (int i = 0; i < N_PARAM; i++)
		{
			uR[i] = PARAMS_REAL[i];
			uI[i] = PARAMS_IMAGINARY[i];
		}
		phiR = PARAM_PHIR;
		phiI = PARAM_PHII;
		//PrintParameters(uR);
		//PrintParameters(uI);
		WriteDataToFile(uR, N_PARAM, "parametersR0", "parameterR");
		WriteDataToFile(uI, N_PARAM, "parametersI0", "parameterI");
	}
	if (processRank == rootRank)
	{
		cout << "PARAMS" << endl;
		PrintParameters(uR);
	}

	sys->InitSystem();
	PostSystemInit();

	int step = 0;
	for (currentTime = 0; currentTime <= TOTALTIME; currentTime += TIMESTEP)
	{
		nAcceptances = 0;
		nTrials = 0;
		step++;
		sys->SetTime(currentTime);

		BroadcastNewParameters(uR, uI, &phiR, &phiI);
		if (USE_MEAN_FOR_FINAL_PARAMETERS)
		{
			if (processRank == rootRank)
			{
				uRList.push_back(ArrayToVector(uR, N_PARAM));
				uIList.push_back(ArrayToVector(uI, N_PARAM));
				phiRList.push_back(phiR);
				phiIList.push_back(phiI);
				if (uRList.size() > 105) //INFO: always keep the last few parameters in a list. at the end of the simulation the average of the last steps is calculated for a final value of the parameter
				{
					uRList.erase(uRList.begin());
					uIList.erase(uIList.begin());
					phiRList.erase(phiRList.begin());
					phiIList.erase(phiIList.begin());
				}
			}
		}
		//PrintParameters(uR);
		//PrintParameters(uI);
		//cout << "phiR=" << phiR << ", phiI=" << phiI << endl;

		AlignCoordinates(R);
		ParallelUpdateExpectationValues(R, uR, uI, phiR, phiI);
		if (processRank == rootRank)
		{
			if (step % 100 == 0)
			{
				WriteDataToFile(localOperators, "localOperators" + to_string(step), "localOperators");
				WriteDataToFile(localEnergyR, "localEnergyR" + to_string(step), "localEnergyR");
				WriteDataToFile(localEnergyI, "localEnergyI" + to_string(step), "localEnergyI");
				WriteDataToFile(localOperatorsMatrix, "localOperatorsMatrix" + to_string(step), "localOperatorsMatrix");
				WriteDataToFile(localOperatorlocalEnergyR, "localOperatorlocalEnergyR" + to_string(step), "localOperatorlocalEnergyR");
				WriteDataToFile(localOperatorlocalEnergyI, "localOperatorlocalEnergyI" + to_string(step), "localOperatorlocalEnergyI");
				WriteDataToFile(otherExpectationValues, "otherExpectationValues" + to_string(step), "Ekin, Epot, wf, g(r)_1, ..., g(r)_100");
			}
			cout << "t=" << currentTime << endl;
			//cout << "localEnergyR=" << localEnergyR << " (" << otherExpectationValues[0] << " + " << otherExpectationValues[1] << ")" << endl;
			cout << "localEnergyR/N=" << localEnergyR / (double)N << " (" << otherExpectationValues[0] / (double)N << " + " << otherExpectationValues[1] / (double)N << ")" << endl;
			AllLocalEnergyR.push_back(localEnergyR);
			AllOtherExpectationValues.push_back(otherExpectationValues);
			AllParametersR.push_back(ArrayToVector(uR, N_PARAM));
			AllParametersI.push_back(ArrayToVector(uI, N_PARAM));
		}
		if (sys->USE_NORMALIZATION_AND_PHASE())
		{
			CalculateNextParameters(R, uR, uI, &phiR, &phiI);
			if (processRank == rootRank)
			{
				//NormalizeWavefunction(otherExpectationValues[2], &phiR);
			}
		}
		else
		{
			CalculateNextParameters(R, uR, uI);
		}
		if (processRank == rootRank)
		{
			if (step % 100 == 0)
			{
				WriteDataToFile(uR, N_PARAM, "parametersR" + to_string(step), "parameterR, phiR=" + to_string(phiR) + ", wf=" + to_string(sys->GetWf()));
				WriteDataToFile(uI, N_PARAM, "parametersI" + to_string(step), "parameterI, phiI=" + to_string(phiI));
			}
			if (step % 100 == 0)
			{
				WriteDataToFile(AllLocalEnergyR, "AllLocalEnergyR", "ER", 100);
				WriteDataToFile(AllOtherExpectationValues, "AllOtherExpectationValues", "kinetic, potential, wf, g(r)", 100);
				WriteDataToFile(AllParametersR, "AllParametersR", "uR", 100);
				WriteDataToFile(AllParametersI, "AllParametersI", "uI", 100);
			}
		}
		if (FileExist("./stop"))
		{
	        cout << "Detected stop-file!" << endl;
	        cout << "Finishing simulation at t=" << currentTime << endl;
	        TOTALTIME = currentTime;
		}
	}

	if (processRank == rootRank)
	{
		PrintParameters(uR);
		PrintParameters(uI);
	}
	//cout << "phiR=" << phiR << ", phiI=" << phiI << endl;

	if (processRank == rootRank)
	{
		Log("Write last files ...");
		WriteDataToFile(AllLocalEnergyR, "AllLocalEnergyR", "ER");
		WriteDataToFile(AllOtherExpectationValues, "AllOtherExpectationValues", "kinetic, potential, wf, g(r)");
		WriteDataToFile(AllParametersR, "AllParametersR", "uR");
		WriteDataToFile(AllParametersI, "AllParametersI", "uI");

		WriteDataToFile(AllLocalEnergyR, "AllLocalEnergyR100", "ER", 100);
		WriteDataToFile(AllOtherExpectationValues, "AllOtherExpectationValues100", "kinetic, potential, wf, g(r)", 100);
		WriteDataToFile(AllParametersR, "AllParametersR100", "uR", 100);
		WriteDataToFile(AllParametersI, "AllParametersI100", "uI", 100);
		WriteRandomGeneratorStatesToFile("state");
		cout << "uniform: " << random01() << endl << "normal: " << randomNormal() << endl;
		cout << "uniform: " << random01() << endl << "normal: " << randomNormal() << endl;
		cout << "uniform: " << random01() << endl << "normal: " << randomNormal() << endl;
	}

	nAcceptances = 0;
	nTrials = 0;
	if (USE_MEAN_FOR_FINAL_PARAMETERS)
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

	AlignCoordinates(R);
	if (processRank == rootRank)
	{
		WriteParticleInputFile("AAFinish_particleconfiguration_" + to_string(N), R);
		cout << "phiR=" << phiR << endl;
		cout << "log(additionalSystemProperties[2])=" << log(additionalSystemProperties[2]) << endl;
		phiR = phiR - log(additionalSystemProperties[2]);
		cout << "phiR=" << phiR << endl;
		WriteConfig("AAFinish_tVMC.config", uR, uI, phiR, phiI);

		if (FileExist("./stop"))
		{
	        cout << "Delete stop-file!" << endl;
	        RemoveFile("./stop");
		}
	}

	sleep(10);
	Log("free memory ...");

	for (int i = 0; i < N; i++)
	{
		delete[] R[i];
	}
	delete[] R;
	delete[] uR;
	delete[] uI;
	delete sys;

	Log("finalize ...");
	MPI_Finalize();

	Log("done");

	return 0;
}

int main(int argc, char **argv) {
	int val = -1;
	cout << "#################################################" << endl;
	cout << "#################################################" << endl;
	Log("start");

	if (argc == 0) //INFO: started with pbs on mach
	{
		string configFilePath = "/home/k3501/k354522/tVMC/bin/tVMC.config";
		val = mainMPI(argc, argv, configFilePath);
	}
	if (argc == 1) //INFO: started without specifying config-file. used for local execution
	{
		string configFilePath = "/home/gartner/Sources/TDVMC/config/drop_3.config";
		//string configFilePath = "/home/gartner/Sources/TDVMC/config/drop_6.config";
		val = mainMPI(argc, argv, configFilePath);
	}
	else if (argc > 1)
	{
		if (strcmp(argv[1], "-test") == 0)
		{
			CalculatePi(argc, argv);
			TestVectorDisplacements();
			val = 0;
		}
		else if (strcmp(argv[1], "-pitest") == 0)
		{
			CalculatePi(argc, argv);
			val = 0;
		}
		else if (strcmp(argv[1], "-vectortest") == 0)
		{
			TestVectorDisplacements();
			val = 0;
		}
		else
		{
			val = mainMPI(argc, argv, string(argv[1]));
		}
	}

	Log("exit");
	return val;
}

