#pragma once

#include "Constants.h"
#include "BulkQT.h"
#include "MathOperators.h"
#include "Utils.h"

#include <iostream>

using namespace std;

class BulkQTPhi: public BulkQT
{

public:
	BulkQTPhi(string configDirectory);

	void CalculateWavefunction(vector<vector<double> >& R, vector<double>& uR, vector<double>& uI, double phiR, double phiI) override;

	void CalculateWFChange(vector<vector<double> >& R, vector<double>& uR, vector<double>& uI, double phiR, double phiI, int changedParticleIndex, vector<double>& oldPosition) override;
};
