#pragma once

#include "Constants.h"
#include "Utils.h"

#include <vector>

using namespace std;

class OneParticleData
{
public:
	vector<double> vecKineticSumR1;
	vector<double> vecKineticSumI1;
	double kineticSumR1;
	double kineticSumI1;
	double kineticSumR1I1;
	double kineticSumR2;
	double kineticSumI2;

public:
	OneParticleData();
};
