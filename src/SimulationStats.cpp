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
	ss << left << setw(wN1) << setfill(c) << "t";
	ss << sep;
	ss << left << setw(wN2) << setfill(c) << "E/N";
	ss << left << setw(wN2) << setfill(c) << "E_kin/N";
	ss << left << setw(wN2) << setfill(c) << "E_pot/N";
	ss << sep;
	ss << left << setw(wN2) << setfill(c) << "exponent";
	ss << left << setw(wN2) << setfill(c) << "wf";
	ss << left << setw(wN2) << setfill(c) << "phiR";
	ss << sep;
	ss << left << setw(wN2) << setfill(c) << "ACC %";
	ss << left << setw(wN2) << setfill(c) << "ACC N";
	ss << sep;
	ss << left << setw(wN2) << setfill(c) << "T/ms (DGL)";
	ss << left << setw(wN2) << setfill(c) << "T/ms (full)";
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
	ss << left << setw(wN2) << setfill(c) << this->locEPot;
	ss << sep;
	ss << left << setw(wN2) << setfill(c) << this->exponent;
	ss << left << setw(wN2) << setfill(c) << this->wf;
	ss << left << setw(wN2) << setfill(c) << this->phiR;
	ss << sep;
	ss << left << setw(wN2) << setfill(c) << this->acceptanceAVGPercent;
	ss << left << setw(wN2) << setfill(c) << this->acceptanceAVGCount;
	//ss << left << setw(wN2) << setfill(c) << this->nTrials;
	ss << sep;
	ss << left << setw(wN2) << setfill(c) << this->durationDGL;
	ss << left << setw(wN2) << setfill(c) << this->durationFullStep;
	string s = ss.str();
	return s;
}

