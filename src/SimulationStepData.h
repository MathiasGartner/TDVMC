#pragma once

#include <vector>

using namespace std;

class SimulationStepData
{
public:
	double time;
	int step;

	double phiR;
	double phiI;
	vector<double> paramsR;
	vector<double> paramsI;

	//difference to parameters in previous timestep eg. param[now] - param[before]
	double phiRDiff;
	double phiIDiff;
	vector<double> paramsRDiff;
	vector<double> paramsIDiff;

	vector<double> localOperators; // for <O_k>
	double localEnergyR; // for <E^R>
	double localEnergyI; // for <E^I>
	vector<vector<double> > localOperatorsMatrix; // for <O_k O_j>
	vector<double> localOperatorlocalEnergyR; // for <O_k E^R>
	vector<double> localOperatorlocalEnergyI; // for <O_k E^I>
	vector<double> otherExpectationValues; // eg. for potential and kinetic energy
};
