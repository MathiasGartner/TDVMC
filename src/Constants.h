#pragma once

#include <math.h>
#include <string>
#include <vector>

using namespace std;

const double HBAR2_2M = 6.0612686;  //hbar^2/2m

extern string OUTPUT_DIRECTORY;
extern int N;           	    	//number of particles
extern int DIM;     	         	//number of dimensions
extern int N_PARAM;	        	  	//number of parameters of trial function
extern double MC_STEP;
extern int MC_NSTEPS;
extern int MC_NTHERMSTEPS;
extern int MC_NINITIALIZATIONSTEPS;
extern int MC_NADDITIONALSTEPS;
extern int MC_NADDITIONALTHERMSTEPS;
extern int MC_NADDITIONALINITIALIZATIONSTEPS;
extern double RHO;
extern double RC;          			//cutoff for WF and LJ
extern double TIMESTEP;
extern double TOTALTIME;
extern int IMAGINARY_TIME;

extern double LBOX;				    //length of box

extern vector<vector<double> > PARAM_SPACE;
