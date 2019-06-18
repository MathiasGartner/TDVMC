#pragma once

#include "BosonClusterWithLog.h"

#include "../Constants.h"
#include "../MathOperators.h"
#include "../Utils.h"

#include <iostream>

using namespace std;

namespace PhysicalSystems
{

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

}
