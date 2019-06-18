#pragma once

#include "BulkQT.h"

#include "../Constants.h"
#include "../MathOperators.h"
#include "../Utils.h"

#include <iostream>

using namespace std;

namespace PhysicalSystems
{

class BulkQTPhi: public BulkQT
{

public:
	BulkQTPhi(vector<double>& params, string configDirectory);

	void CalculateWavefunction(vector<vector<double> >& R, vector<double>& uR, vector<double>& uI, double phiR, double phiI) override;

	void CalculateWFChange(vector<vector<double> >& R, vector<double>& uR, vector<double>& uI, double phiR, double phiI, int changedParticleIndex, vector<double>& oldPosition) override;
};

}
