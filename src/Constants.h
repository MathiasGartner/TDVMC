/*-----------------------------------------------------------------------------
 *
 * 		Name:			Constants.h
 * 		Author:			Mathias Gartner
 * 		Description:	Constants for the simulation that are shared within the
 * 						whole tVMC program.
 *
 * --------------------------------------------------------------------------*/

#pragma once

#include <math.h>
#include <string>
#include <vector>

using namespace std;

//const double HBAR2_2M = 6.0612686;  //hbar^2/2m for 4He droplets
//const double HBAR2_2M = 6.06359;  //hbar^2/2m for X_4 cluster from split
//const double HBAR2_2M = 0.5;  //hbar^2/2m for model systems
const double HBAR2_2M = 1.0;  //hbar^2/2m for model systems

extern string OUT_DIR;
extern int N;           	    	//number of particles
extern int DIM;     	         	//number of dimensions
extern int N_PARAM;	        	  	//number of parameters of trial function
extern int IMAGINARY_TIME;					//0: real time, 1: imaginary time

extern double LBOX;				    //length of box
extern double LBOX_R;				//reciprocal length of box
extern double LBOX_2;				//half of box length

namespace Constants
{
	namespace Split
	{
		const double hbar = 1.054571628e-34;
		const double u = 1.660538782e-27;
		const double mK2K = 0.001;
		const double A2m = 1e-10;
		const double kb = 1.3806504e-23;
	}
}
