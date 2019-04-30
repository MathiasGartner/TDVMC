#pragma once

#include <math.h>
#include <string>
#include <vector>

using namespace std;

//const double HBAR2_2M = 6.0612686;  //hbar^2/2m for 4He droplets
const double HBAR2_2M = 6.06359;  //hbar^2/2m for X_4 cluster from split

extern string OUT_DIR;
extern int N;           	    	//number of particles
extern int DIM;     	         	//number of dimensions
extern int N_PARAM;	        	  	//number of parameters of trial function

extern double LBOX;				    //length of box
extern double LBOX_R;				//reciprocal length of box
