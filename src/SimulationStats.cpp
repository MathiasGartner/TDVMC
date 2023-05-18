/*-----------------------------------------------------------------------------
 *
 * 		Name:			SimulationStats.cpp
 * 		Author:			Mathias Gartner
 * 		Description:	see SimulationStats.h
 *
 * --------------------------------------------------------------------------*/

#include "SimulationStats.h"

string SimulationStats::getHeader()
{
	stringstream ss;
	//groups
	ss << left << setw(wN1) << setfill(c) << "## time";
	ss << sep;
	ss << left << setw(wN2 + wN2 + wN2 + wN2) << setfill(c) << "## energy";
	ss << sep;
	ss << left << setw(wN2 + wN2 + wN2) << setfill(c) << "## wavefunction";
	ss << sep;
	ss << left << setw(wN2 + wN2) << setfill(c) << "## sampling acceptance";
	ss << sep;
	ss << left << setw(wN2 + wN2) << setfill(c) << "## runtimes in ms";
	ss << endl;
	//fields
	ss << left << setw(wN1) << setfill(c) << "t";
	ss << sep;
	ss << left << setw(wN2) << setfill(c) << "E/N";
	ss << left << setw(wN2) << setfill(c) << "E_kin/N";
	ss << left << setw(wN2) << setfill(c) << "E_ext/N";
	ss << left << setw(wN2) << setfill(c) << "E_int/N";
	ss << sep;
	ss << left << setw(wN2) << setfill(c) << "exponent";
	ss << left << setw(wN2) << setfill(c) << "wf";
	ss << left << setw(wN2) << setfill(c) << "phiR";
	ss << sep;
	ss << left << setw(wN2) << setfill(c) << "ACC %";
	ss << left << setw(wN2) << setfill(c) << "ACC N";
	ss << sep;
	ss << left << setw(wN1) << setfill(c) << "MC";
	ss << left << setw(wN1) << setfill(c) << "EOM";
	ss << left << setw(wN1) << setfill(c) << "full step";
	string s = ss.str();
	return s;
}

string SimulationStats::toString()
{
	stringstream ss;
	ss << left << setw(wN1) << setfill(c) << this->currentTime;
	ss << sep;
	ss << left << setw(wN2) << setfill(c) << this->locE;
	ss << left << setw(wN2) << setfill(c) << this->locEKin;
	ss << left << setw(wN2) << setfill(c) << this->locEExt;
	ss << left << setw(wN2) << setfill(c) << this->locEInt;
	ss << sep;
	ss << left << setw(wN2) << setfill(c) << this->exponent;
	ss << left << setw(wN2) << setfill(c) << this->wf;
	ss << left << setw(wN2) << setfill(c) << this->phiR;
	ss << sep;
	ss << left << setw(wN2) << setfill(c) << this->acceptanceAVGPercent;
	ss << left << setw(wN2) << setfill(c) << this->acceptanceAVGCount;
	//ss << left << setw(wN2) << setfill(c) << this->nTrials;
	ss << sep;
	ss << left << setw(wN1) << setfill(c) << this->durationSampling;
	ss << left << setw(wN1) << setfill(c) << this->durationSolveEOM;
	ss << left << setw(wN1) << setfill(c) << this->durationFullStep;
	string s = ss.str();
	return s;
}

