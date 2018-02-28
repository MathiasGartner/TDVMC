/*
 ============================================================================
 Name        : TDVMC.c
 Author      : Mathias Gartner
 Version     :
 Copyright   : 
 Description :
 ============================================================================
 */

#include <chrono>
#include "Constants.h"
#include <cstring>
#include <ctime>
#include <fstream>
#include "HeCubic5.h"
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

//bool useNIC = false; //for drops
bool useNIC = true; //for bulk with PBC

string OUTPUT_DIRECTORY;
string originalOutputDirectory;
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
vector<vector<double> > AllAdditionalSystemProperties;

default_random_engine generator;
uniform_real_distribution<double> distUniform(0.0, 1.0);
normal_distribution<double> distNormal(0.0,1.0);

int main(int argc, char **argv) {
	int val = -1;
	cout << "start" << endl;

	if (argc > 1)
	{
		if (strcmp(argv[1], "-pitest") == 0)
		{
			CalculatePi(argc, argv);
		}
		else
		{
			//TODO: TDVMC
		}
	}

	val = 0;
	cout << "finish " << val << endl;
	return val;
}

