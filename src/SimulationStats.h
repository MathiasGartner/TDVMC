/*-----------------------------------------------------------------------------
 *
 * 		Name:			SimulationStats.h
 * 		Author:			Mathias Gartner
 * 		Description:	Holds information about the current status of the
 * 						simulation. Used for printing status reports.
 *
 * --------------------------------------------------------------------------*/

#pragma once

#include <iomanip>
#include <sstream>
#include <string>

using namespace std;

class SimulationStats
{
public:
	long long nTrials;
	long long acceptanceAVGCount;
	double acceptanceAVGPercent;

	double currentTime;

	double locE;
	double locEKin;
	double locEPot;

	double exponent;
	double wf;
	double phiR;

	double durationDGL;
	double durationFullStep;


	int wN1 = 9;
	int wN2 = 12;
	int wN3 = 14;
	char c = ' ';
	string sep = "| ";

	string getHeader();
	string toString();
};
