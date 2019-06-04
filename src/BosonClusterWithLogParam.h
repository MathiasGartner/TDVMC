#pragma once

#include "Constants.h"
#include "ICorrelatedSamplingData.h"
#include "BosonClusterWithLog.h"
#include "MathOperators.h"
#include "Utils.h"

#include <iostream>

using namespace std;

class BosonClusterWithLogParam: public BosonClusterWithLog
{
private:

public:
	BosonClusterWithLogParam(vector<double>& params, string configDirectory);

	void InitSystem();

	void CalculateExpectationValues(vector<vector<double> >& R, vector<double>& uR, vector<double>& uI, double phiR, double phiI);
	void CalculateWavefunction(vector<vector<double> >& R, vector<double>& uR, vector<double>& uI, double phiR, double phiI);
	void CalculateWFChange(vector<vector<double> >& R, vector<double>& uR, vector<double>& uI, double phiR, double phiI, int changedParticleIndex, vector<double>& oldPosition);
};

