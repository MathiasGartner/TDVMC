#include "mpi.h"

#include "BulkSplines.h"
#include "BulkSplinesPhi.h"
#include "BulkQT.h"
#include "BulkQTPhi.h"
#include "ConfigItem.h"
#include "Constants.h"
#include "GaussianWavepacket.h"
#include "HeBulk.h"
#include "HeDrop.h"
#include "MathOperators.h"
#include "MPIMethods.h"
#include "Timer.h"
#include "Utils.h"
#include "VMCSampler.h"

#include "test/MPITest.h"
#include "test/PiCalculator.h"
#include "test/Tests.h"

//#undef SEEK_SET
//#undef SEEK_END
//#undef SEEK_CUR

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

#include <json/json.h>

using namespace std;

IPhysicalSystem* sys;

vector<ConfigItem> configItems;

string SYSTEM_TYPE;
string OUTPUT_DIRECTORY;		//from config file
string configDirectory;
string configFilePath;
string OUT_DIR;					//directory name generated from parameter settings
int N;           	    		//number of particles
int DIM;     	        	 	//number of dimensions
int N_PARAM;	        	  	//number of parameters of trial function
double MC_STEP;
double MC_STEP_OFFSET;
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
int ODE_SOLVER_TYPE;
double LBOX;
double LBOX_R;
vector<double> PARAMS_REAL;
vector<double> PARAMS_IMAGINARY;
double PARAM_PHIR;
double PARAM_PHII;
int USE_PARAMETER_ACCEPTANCE_CHECK;
int PARAMETER_ACCEPTANCE_CHECK_TYPE;
int WRITE_EVERY_NTH_STEP_TO_FILE;
int MC_NSTEP_MULTIPLICATION_FACTOR_FOR_WRITE_DATA;
int USE_MEAN_FOR_FINAL_PARAMETERS;
int USE_NORMALIZE_WF;
int USE_ADJUST_PARAMETERS;

int numOfProcesses = 1;
int rootRank = 0;
int processRank = 0;
int nAcceptances = 0;
int nTrials = 0;
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

vector<vector<double> > uRList;
vector<vector<double> > uIList;
vector<double> phiRList;
vector<double> phiIList;
vector<vector<double> > uRListDiffs;
vector<vector<double> > uIListDiffs;
vector<double> phiRListDiffs;
vector<double> phiIListDiffs;

vector<double> AllLocalEnergyR;
vector<vector<double> > AllOtherExpectationValues;
vector<vector<double> > AllParametersR;
vector<vector<double> > AllParametersI;
vector<vector<double> > AllAdditionalSystemProperties;

vector<vector<vector<double> > > mcSamples;

mt19937_64 generator;
//default_random_engine generator;
uniform_real_distribution<double> distUniform(0.0, 1.0);
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

void PrintData(vector<double> data)
{
	for (unsigned int i = 0; i < data.size(); i++)
	{
		cout << data[i] << ", ";
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

double randomNormal(double sigma, double mu)
{
	return randomNormal() * sigma + mu;
}

int randomInt(int maxValue)
{
	//TODO: test for maxValue > 0 is only needed for Gaussian wavepacket simulation where only one  particle is used
	return maxValue > 0 ? ((int) floor(random01() * maxValue)) % maxValue : maxValue;
}

int randomInt(int minValue, int maxValue)
{
	return randomInt(maxValue - minValue) + minValue;
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
}

///////////////////////////
/// configuration items ///
///////////////////////////

void RegisterAllConfigItems()
{
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
	configItems.push_back(ConfigItem("MC_NADDITIONALSTEPS", &MC_NADDITIONALSTEPS, ConfigItemType::INT));
	configItems.push_back(ConfigItem("MC_NADDITIONALTHERMSTEPS", &MC_NADDITIONALTHERMSTEPS, ConfigItemType::INT));
	configItems.push_back(ConfigItem("MC_NADDITIONALINITIALIZATIONSTEPS", &MC_NADDITIONALINITIALIZATIONSTEPS, ConfigItemType::INT));
	configItems.push_back(ConfigItem("TIMESTEP", &TIMESTEP, ConfigItemType::DOUBLE));
	configItems.push_back(ConfigItem("TOTALTIME", &TOTALTIME, ConfigItemType::DOUBLE));
	configItems.push_back(ConfigItem("IMAGINARY_TIME", &IMAGINARY_TIME, ConfigItemType::INT));
	configItems.push_back(ConfigItem("ODE_SOLVER_TYPE", &ODE_SOLVER_TYPE, ConfigItemType::INT));
	configItems.push_back(ConfigItem("USE_PARAMETER_ACCEPTANCE_CHECK", &USE_PARAMETER_ACCEPTANCE_CHECK, ConfigItemType::INT));
	configItems.push_back(ConfigItem("PARAMETER_ACCEPTANCE_CHECK_TYPE", &PARAMETER_ACCEPTANCE_CHECK_TYPE, ConfigItemType::INT));
	configItems.push_back(ConfigItem("WRITE_EVERY_NTH_STEP_TO_FILE", &WRITE_EVERY_NTH_STEP_TO_FILE, ConfigItemType::INT));
	configItems.push_back(ConfigItem("MC_NSTEP_MULTIPLICATION_FACTOR_FOR_WRITE_DATA", &MC_NSTEP_MULTIPLICATION_FACTOR_FOR_WRITE_DATA, ConfigItemType::INT));
	configItems.push_back(ConfigItem("PARAMS_REAL", PARAMS_REAL, ConfigItemType::ARR_DOUBLE));
	configItems.push_back(ConfigItem("PARAMS_IMAGINARY", PARAMS_IMAGINARY, ConfigItemType::ARR_DOUBLE));
	configItems.push_back(ConfigItem("PARAM_PHIR", &PARAM_PHIR, ConfigItemType::DOUBLE));
	configItems.push_back(ConfigItem("PARAM_PHII", &PARAM_PHII, ConfigItemType::DOUBLE));
	configItems.push_back(ConfigItem("USE_MEAN_FOR_FINAL_PARAMETERS", &USE_MEAN_FOR_FINAL_PARAMETERS, ConfigItemType::INT));
	configItems.push_back(ConfigItem("USE_NORMALIZE_WF", &USE_NORMALIZE_WF, ConfigItemType::INT));
	configItems.push_back(ConfigItem("USE_ADJUST_PARAMETERS", &USE_ADJUST_PARAMETERS, ConfigItemType::INT));
}

void ReadConfig(string filePath)
{
	//cout << "file exists: " << FileExist(filePath) << "." << endl;
	Json::Value configData;
	Json::Reader configReader;
	ifstream configFile(filePath, ifstream::binary);
	configReader.parse(configFile, configData, false);

	for (auto ci : configItems)
	{
		ci.setValue(configData[ci.name]);
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
	for (auto ci : configItems)
	{
		configFile << "\t" << "\"" << ci.name << "\"" << " : " << ci.getJsonString() << "," << endl;
	}
	configFile << "}";

	configFile.close();
}

void BroadcastConfig()
{
	for (auto ci : configItems)
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
		srand(processRank + 1);
	}
	LBOX = pow((N / RHO), 1.0 / DIM); //box with dimensions [-L/2, L/2]
	LBOX_R = 1.0 / LBOX;
	nAcceptances = 0;
	nTrials = 0;

	mc_nsteps = (double) MC_NSTEPS;
	mc_nsteps_original = MC_NSTEPS;
	mc_nadditionalsteps = (double) MC_NADDITIONALSTEPS;

	mcSamples.resize(MC_NSTEPS);
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
	return successful;
}

void InitCoordinateConfiguration(vector<vector<double> >& R)
{
	int fileFound = 0;
	if (processRank == rootRank)
	{
		string filename = "particleconfiguration_" + to_string(N);
		if (LoadLastPositionsFromFile(filename, R))
		{
			fileFound = 1;
		}
	}
	MPIMethods::BroadcastValue(&fileFound);
	if (fileFound)
	{
		MPIMethods::BroadcastValues(R);
	}
	else
	{
		if (processRank == rootRank)
		{
			cout << "LBOX=" << LBOX << endl;
		}
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
			double n = round(pow(N, 1.0 / 3.0));
			double l = pow(LBOX, 1.0 / 3.0);
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
	double p;
	bool sampleOkay = true;
	int randomParticle = randomInt(N - 1);
	vector<double> oldPosition(R[randomParticle]);
	for (int i = 0; i < DIM; i++)
	{
		R[randomParticle][i] += randomNormal(MC_STEP, MC_STEP_OFFSET);
	}
	double wfQuotient = sys->CalculateWFQuotient(R, uR, uI, phiR, phiI, randomParticle, oldPosition);

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
			mcSamples[i] = R;
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

void UpdateExpectationValuesForGivenSamples(vector<vector<vector<double> > >& samples, vector<double>& uR, vector<double>& uI, double phiR, double phiI)
{
	int count = samples.size();
	double dblCount = (double) count;

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

		//INFO: calculate contribution to average value
		localOperators += sys->GetLocalOperators() / dblCount;
		localEnergyR += sys->GetLocalEnergyR() / dblCount;
		localEnergyI += sys->GetLocalEnergyI() / dblCount;
		localOperatorsMatrix += sys->GetLocalOperatorsMatrix() / dblCount;
		localOperatorlocalEnergyR += sys->GetLocalOperatorlocalEnergyR() / dblCount;
		localOperatorlocalEnergyI += sys->GetLocalOperatorlocalEnergyI() / dblCount;
		otherExpectationValues += sys->GetOtherExpectationValues() / dblCount;

		if ((10 * i) % count == 0)
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
}

void ParallelUpdateExpectationValuesForGivenSamples(vector<vector<vector<double> > >& samples, vector<double>& uR, vector<double>& uI, double phiR, double phiI)
{
	//Timer t;
	//double dblDuration;
	//t.start();

	UpdateExpectationValuesForGivenSamples(samples, uR, uI, phiR, phiI);

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

void NormalizeParameters(vector<double>& uR, vector<double>& uI, double *phiR, double *phiI)
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
	if (sys->USE_NORMALIZATION_AND_PHASE)
	{
		BuildSystemOfEquationsForParametersIncludePhi(matrix, energiesReal, energiesImag);
	}
	else
	{
		//BuildSystemOfEquationsForParametersNoPhi(matrix, energiesReal, energiesImag);
		BuildSystemOfEquationsForParametersIncludePhi(matrix, energiesReal, energiesImag);
	}
}

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
					logged = true;
				}
			}
		}
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

void CalculateNextParametersEuler(vector<double>& uR, vector<double>& uI, double *phiR, double *phiI)
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

void CalculateNextParametersPC(vector<vector<double> >& R, vector<double>& uR, vector<double>& uI, double *phiR, double *phiI)
{
	vector<double> energiesReal; //rhs of the equation that corresponds to uDotR;
	vector<double> energiesImag; //rhs of the equation that corresponds to uDotI;
	vector<vector<double> > matrix;
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

	//TODO: is this initialization needed?
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
}

void CalculateNextParametersRK4(vector<vector<double> >& R, vector<double>& uR, vector<double>& uI, double *phiR, double *phiI)
{
	vector<double> energiesReal; //rhs of the equation that corresponds to uDotR;
	vector<double> energiesImag; //rhs of the equation that corresponds to uDotI;
	vector<vector<double> > matrix;

	vector<double> tmpUR(N_PARAM);
	vector<double> tmpUI(N_PARAM);
	double tmpPhiR = 0;
	double tmpPhiI = 0;

	vector<vector<double> > uDotR;
	vector<vector<double> > uDotI;
	vector<double> phiDotR;
	vector<double> phiDotI;

	//TODO: is this initialization needed?
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
			tmpUI[i] = uI[i] + uDotI[0][i] * TIMESTEP / 2.0;
		}
		tmpPhiR = *phiR + phiDotR[0] * TIMESTEP / 2.0;
		tmpPhiI = *phiI + phiDotI[0] * TIMESTEP / 2.0;
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
			tmpUI[i] = uI[i] + uDotI[1][i] * TIMESTEP / 2.0;
		}
		tmpPhiR = *phiR + phiDotR[1] * TIMESTEP / 2.0;
		tmpPhiI = *phiI + phiDotI[1] * TIMESTEP / 2.0;
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
			tmpUI[i] = uI[i] + uDotI[2][i] * TIMESTEP;
		}
		tmpPhiR = *phiR + phiDotR[2] * TIMESTEP;
		tmpPhiI = *phiI + phiDotI[2] * TIMESTEP;
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
			uI[i] = uI[i] + ((uDotI[0][i] + 2.0 * uDotI[1][i] + 2.0 * uDotI[2][i] + uDotI[3][i]) / 6.0) * TIMESTEP;
		}
		*phiR = *phiR + ((phiDotR[0] + 2.0 * phiDotR[1] + 2.0 * phiDotR[2] + phiDotR[3]) / 6.0) * TIMESTEP;
		*phiI = *phiI + ((phiDotI[0] + 2.0 * phiDotI[1] + 2.0 * phiDotI[2] + phiDotI[3]) / 6.0) * TIMESTEP;
	}
}

void CalculateNextParametersRK4ReuseSamples(vector<vector<double> >& R, vector<double>& uR, vector<double>& uI, double *phiR, double *phiI)
{
	vector<double> energiesReal; //rhs of the equation that corresponds to uDotR;
	vector<double> energiesImag; //rhs of the equation that corresponds to uDotI;
	vector<vector<double> > matrix;

	vector<double> tmpUR(N_PARAM);
	vector<double> tmpUI(N_PARAM);
	double tmpPhiR = 0;
	double tmpPhiI = 0;

	vector<vector<double> > uDotR;
	vector<vector<double> > uDotI;
	vector<double> phiDotR;
	vector<double> phiDotI;

	//TODO: is this initialization needed?
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
			tmpUI[i] = uI[i] + uDotI[0][i] * TIMESTEP / 2.0;
		}
		tmpPhiR = *phiR + phiDotR[0] * TIMESTEP / 2.0;
		tmpPhiI = *phiI + phiDotI[0] * TIMESTEP / 2.0;
	}
	BroadcastNewParameters(tmpUR, tmpUI, &tmpPhiR, &tmpPhiI);
	ParallelUpdateExpectationValuesForGivenSamples(mcSamples, tmpUR, tmpUI, tmpPhiR, tmpPhiI);

	if (processRank == rootRank)
	{
		BuildSystemOfEquationsForParameters(matrix, energiesReal, energiesImag);
		PerformCholeskyDecomposition(matrix);
		SolveForParametersDot(matrix, energiesReal, energiesImag, uDotR[1], uDotI[1], &(phiDotR[1]), &(phiDotI[1]));
		for (int i = 0; i < N_PARAM; i++)
		{
			tmpUR[i] = uR[i] + uDotR[1][i] * TIMESTEP / 2.0;
			tmpUI[i] = uI[i] + uDotI[1][i] * TIMESTEP / 2.0;
		}
		tmpPhiR = *phiR + phiDotR[1] * TIMESTEP / 2.0;
		tmpPhiI = *phiI + phiDotI[1] * TIMESTEP / 2.0;
	}
	BroadcastNewParameters(tmpUR, tmpUI, &tmpPhiR, &tmpPhiI);
	ParallelUpdateExpectationValuesForGivenSamples(mcSamples, tmpUR, tmpUI, tmpPhiR, tmpPhiI);

	if (processRank == rootRank)
	{
		BuildSystemOfEquationsForParameters(matrix, energiesReal, energiesImag);
		PerformCholeskyDecomposition(matrix);
		SolveForParametersDot(matrix, energiesReal, energiesImag, uDotR[2], uDotI[2], &(phiDotR[2]), &(phiDotI[2]));
		for (int i = 0; i < N_PARAM; i++)
		{
			tmpUR[i] = uR[i] + uDotR[2][i] * TIMESTEP;
			tmpUI[i] = uI[i] + uDotI[2][i] * TIMESTEP;
		}
		tmpPhiR = *phiR + phiDotR[2] * TIMESTEP;
		tmpPhiI = *phiI + phiDotI[2] * TIMESTEP;
	}
	BroadcastNewParameters(tmpUR, tmpUI, &tmpPhiR, &tmpPhiI);
	ParallelUpdateExpectationValuesForGivenSamples(mcSamples, tmpUR, tmpUI, tmpPhiR, tmpPhiI);

	if (processRank == rootRank)
	{
		BuildSystemOfEquationsForParameters(matrix, energiesReal, energiesImag);
		PerformCholeskyDecomposition(matrix);
		SolveForParametersDot(matrix, energiesReal, energiesImag, uDotR[3], uDotI[3], &(phiDotR[3]), &(phiDotI[3]));
		for (int i = 0; i < N_PARAM; i++)
		{
			uR[i] = uR[i] + ((uDotR[0][i] + 2.0 * uDotR[1][i] + 2.0 * uDotR[2][i] + uDotR[3][i]) / 6.0) * TIMESTEP;
			uI[i] = uI[i] + ((uDotI[0][i] + 2.0 * uDotI[1][i] + 2.0 * uDotI[2][i] + uDotI[3][i]) / 6.0) * TIMESTEP;
		}
		*phiR = *phiR + ((phiDotR[0] + 2.0 * phiDotR[1] + 2.0 * phiDotR[2] + phiDotR[3]) / 6.0) * TIMESTEP;
		*phiI = *phiI + ((phiDotI[0] + 2.0 * phiDotI[1] + 2.0 * phiDotI[2] + phiDotI[3]) / 6.0) * TIMESTEP;
	}
}

void CalculateNextParameters(vector<vector<double> >& R, vector<double>& uR, vector<double>& uI, double *phiR, double *phiI)
{
	Timer t;
	if (processRank == rootRank)
	{
		t.start();
	}

	if (ODE_SOLVER_TYPE == 0)
	{
		CalculateNextParametersEuler(uR, uI, phiR, phiI);
	}
	else if (ODE_SOLVER_TYPE == 1)
	{
		vector<vector<double> >& Rcopy(R);
		CalculateNextParametersPC(Rcopy, uR, uI, phiR, phiI);
	}
	else if (ODE_SOLVER_TYPE == 2)
	{
		vector<vector<double> >& Rcopy(R);
		CalculateNextParametersRK4(Rcopy, uR, uI, phiR, phiI);
	}
	else if (ODE_SOLVER_TYPE == 3)
	{
		vector<vector<double> >& Rcopy(R);
		CalculateNextParametersRK4ReuseSamples(Rcopy, uR, uI, phiR, phiI);
	}

	if (processRank == rootRank)
	{
		t.stop();
		Log("DGL duration = " + to_string(t.duration()) + " ms");
	}
}

void CalculateNextParameters(vector<vector<double> >& R, vector<double>& uR, vector<double>& uI)
{
	double tmpPhiR = 0;
	double tmpPhiI = 0;
	CalculateNextParameters(R, uR, uI, &tmpPhiR, &tmpPhiI);
}

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
bool AcceptNewParams(vector<double>& uR, vector<double>& uI, double phiR, double phiI)
{
	bool accept = true;
	if (PARAMETER_ACCEPTANCE_CHECK_TYPE == 0)
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
		double threshold = 0.05;
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
	return accept;
}

int mainMPI(int argc, char** argv)
{
	char processName[80];
	int processNameLength;

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &processRank);
	MPI_Comm_size(MPI_COMM_WORLD, &numOfProcesses);
	MPI_Get_processor_name(processName, &processNameLength);

	MPIMethods::numOfProcesses = numOfProcesses;
	MPIMethods::processRank = processRank;
	MPIMethods::rootRank = rootRank;

	//Log("running on cpu " + to_string(get_cpu_id()));
	MPIMethods::GetCPUAllocation();

	RegisterAllConfigItems();
	if (processRank == rootRank)
	{
		Log("Master process started...");

		if (FileExist(configFilePath))
		{
			ReadConfig(configFilePath);
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

	if (SYSTEM_TYPE == "HeDrop")
	{
		sys = new HeDrop(configDirectory);
	}
	else if (SYSTEM_TYPE == "HeBulk")
	{
		sys = new HeBulk(configDirectory);
	}
	else if (SYSTEM_TYPE == "BulkSplines")
	{
		//sys = new BulkSplines(configDirectory);
		sys = new BulkSplinesPhi(configDirectory);
	}
	else if (SYSTEM_TYPE == "BulkQT")
	{
		//sys = new BulkQT(configDirectory);
		sys = new BulkQTPhi(configDirectory);
	}
	else if (SYSTEM_TYPE == "GaussianWavepacket")
	{
		sys = new GaussianWavepacket(configDirectory);
	}
	else
	{
		if (processRank == rootRank)
		{
			Log("System type \"" + SYSTEM_TYPE + "\" not available.", ERROR);
		}
		return 1;
	}
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

	for (unsigned int i = 0; i < R.size(); i++)
	{
		R[i].resize(DIM);
	}

	InitCoordinateConfiguration(R);
	if (processRank == rootRank)
	{
		WriteParticleInputFile(OUT_DIR + "particleconfiguration.csv", R);
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

		if (processRank == rootRank)// && USE_ADJUST_PARAMETERS == true)
		{
			//AdjustParameters(uR, uI, &phiR, &phiI);
		}

		WriteDataToFile(uR, "parametersR0", "parameterR");
		WriteDataToFile(uI, "parametersI0", "parameterI");
	}

	sys->InitSystem();
	PostSystemInit();

	int step = 0;
	int acceptNewParams;
	int nrOfAcceptParameterTrials;

	for (currentTime = 0; currentTime <= TOTALTIME; currentTime += TIMESTEP)
	{
		acceptNewParams = 0;
		nrOfAcceptParameterTrials = 0;
		MC_NSTEPS = mc_nsteps_original;
		nAcceptances = 0;
		nTrials = 0;
		step++;
		sys->SetTime(currentTime);
		sys->SetStep(step);

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

		while (acceptNewParams != 1 && nrOfAcceptParameterTrials < 4)
		{
			nrOfAcceptParameterTrials++;
			MC_NSTEPS *= nrOfAcceptParameterTrials;
			AlignCoordinates(R);
			ParallelUpdateExpectationValues(R, uR, uI, phiR, phiI);
			if (processRank == rootRank)
			{
				if (step % WRITE_EVERY_NTH_STEP_TO_FILE == 0)
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
				//cout << "localEnergyR=" << localEnergyR << " (" << otherExpectationValues[0] << " + " << otherExpectationValues[1] << ")" << endl;
				cout << "localEnergyR/N=" << localEnergyR / (double) N << " (" <<
						otherExpectationValues[0] / (double) N << " + " <<
						otherExpectationValues[1] / (double) N << " + " <<
						//otherExpectationValues[2] / (double) N << " + " <<
						//otherExpectationValues[3] / (double) N << ")" <<
						endl;
				//cout << "localEnergyR/N=" << localEnergyR / (double)N << endl;
				AllLocalEnergyR.push_back(localEnergyR);
				AllOtherExpectationValues.push_back(otherExpectationValues);
				AllParametersR.push_back(uR);
				AllParametersR[AllParametersR.size() - 1].push_back(phiR);
				AllParametersI.push_back(uI);
				AllParametersI[AllParametersI.size() - 1].push_back(phiI);
			}
			if (sys->USE_NORMALIZATION_AND_PHASE)
			{
				CalculateNextParameters(R, uR, uI, &phiR, &phiI);
				if (processRank == rootRank && USE_NORMALIZE_WF == 1)
				{
					NormalizeWavefunction(sys->GetWf(), &phiR);
				}
			}
			else
			{
				CalculateNextParameters(R, uR, uI);
			}
			if (processRank == rootRank && USE_ADJUST_PARAMETERS == 1)
			{
				AdjustParameters(uR, uI, &phiR, &phiI);
			}
			if (USE_PARAMETER_ACCEPTANCE_CHECK == 1)
			{
				if (processRank == rootRank)
				{
					acceptNewParams = AcceptNewParams(uR, uI, phiR, phiI) ? 1 : 0;
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

		if (processRank == rootRank)
		{
			if (step % WRITE_EVERY_NTH_STEP_TO_FILE == 0)
			{
				WriteDataToFile(uR, "parametersR" + to_string(step), "parameterR, phiR=" + to_string(phiR) + ", wf=" + to_string(sys->GetWf()));
				WriteDataToFile(uI, "parametersI" + to_string(step), "parameterI, phiI=" + to_string(phiI));
			}
			if (step % WRITE_EVERY_NTH_STEP_TO_FILE == 0)
			{
				WriteDataToFile(AllLocalEnergyR, "AllLocalEnergyR", "ER", WRITE_EVERY_NTH_STEP_TO_FILE);
				WriteDataToFile(AllOtherExpectationValues, "AllOtherExpectationValues", "kinetic, potential, wf, g(r)", WRITE_EVERY_NTH_STEP_TO_FILE);
				WriteDataToFile(AllParametersR, "AllParametersR", "uR", WRITE_EVERY_NTH_STEP_TO_FILE);
				WriteDataToFile(AllParametersI, "AllParametersI", "uI", WRITE_EVERY_NTH_STEP_TO_FILE);
			}
		}

		//INFO: check if simulation should be cancelled
		int cancel = 0;
		if (processRank == rootRank)
		{
			if (FileExist("./stop"))
			{
				Log("Detected stop-file!", WARNING);
				cancel = 1;
			}
			if (!isfinite(localEnergyR))
			{
				Log("Energy not finite", ERROR);
				cancel = 2;
			}
			if (cancel != 0)
			{
				Log("Finishing simulation at t=" + to_string(currentTime));
			}
		}
		MPIMethods::BroadcastValue(&cancel);
		if (cancel != 0)
		{
			TOTALTIME = currentTime;
		}
	}

	if (processRank == rootRank)
	{
		PrintData(uR);
		PrintData(uI);
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
	}
	if (processRank == rootRank)
	{
		int tmp = 0;
		tmp = system(("mkdir " + OUT_DIR + "random").c_str());
		cout << tmp << endl;
	}
	MPIMethods::Barrier();
	WriteRandomGeneratorStatesToFile("random/state");
	if (processRank == rootRank)
	{
		cout << "uniform: " << random01() << endl << "normal: " << randomNormal() << endl;
		cout << "uniform: " << random01() << endl << "normal: " << randomNormal() << endl;
		cout << "uniform: " << random01() << endl << "normal: " << randomNormal() << endl;
	}

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

	AlignCoordinates(R);
	if (processRank == rootRank)
	{
		MC_NSTEPS = mc_nsteps_original;
		WriteParticleInputFile("AAFinish_particleconfiguration_" + to_string(N), R);
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

	Log("free memory ...");
	//TODO: how to properly delete the IPhysicalSystem pointer?
	//delete sys;

	Log("finalize ...");
	MPIMethods::Barrier();
	MPI_Finalize();

	Log("done");

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
		energies.push_back({e, -VMCSampler::energy1 / (double)N, -VMCSampler::energy2 / (double)N, -VMCSampler::energy2_1 / (double)N, -VMCSampler::energy2_2 / (double)N});
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

	if (argc == 0) //INFO: started with pbs on mach
	{
		configFilePath = "/home/k3501/k354522/tVMC/bin/tVMC.config";
		val = mainMPI(argc, argv);
	}
	else if (argc == 1) //INFO: started without specifying config-file. used for local execution
	{
		//SYSTEM_TYPE = "BulkQT";
		SYSTEM_TYPE = "BulkSplines";
		if (SYSTEM_TYPE == "HeDrop")
		{
			configFilePath = "/home/gartner/Sources/TDVMC/config/drop_3.config";
		}
		else if (SYSTEM_TYPE == "HeBulk")
		{
			configFilePath = "/home/gartner/Sources/TDVMC/config/bulk_64.config";
		}
		else if (SYSTEM_TYPE == "BulkSplines")
		{
			configFilePath = "/home/gartner/Sources/TDVMC/config/bulkSplines.config";
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
		if (processRank == rootRank)
		{
			Log("configFilePath: " + configFilePath);
		}
		val = mainMPI(argc, argv);
	}
	else if (argc > 1)
	{
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
		}
	}

	//Log("exit");
	return val;
}

