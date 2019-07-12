#pragma once

#include "BosonCluster.h"

#include "../Constants.h"
#include "../MathOperators.h"
#include "../Utils.h"

#include <iostream>

using namespace std;

namespace PhysicalSystems
{

class BosonClusterWithLog: public BosonCluster
{
private:
	double logPrefactor;

protected:
	double logSum;
	vector<vector<double> > logSumD;
	vector<double> logSumD2;
	double logSumNew;

public:
	BosonClusterWithLog(vector<double>& params, string configDirectory);

	void InitSystem() override;

	void CalculateExpectationValues(vector<vector<double> >& R, vector<double>& uR, vector<double>& uI, double phiR, double phiI) override;

	void CalculateWavefunction(vector<vector<double> >& R, vector<double>& uR, vector<double>& uI, double phiR, double phiI) override;

	void CalculateWFChange(vector<vector<double> >& R, vector<double>& uR, vector<double>& uI, double phiR, double phiI, int changedParticleIndex, vector<double>& oldPosition) override;

	void AcceptMove() override;
};

}
