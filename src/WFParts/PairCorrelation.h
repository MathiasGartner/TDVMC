#pragma once

#include "SplinedFunction.h"

using namespace std;

namespace WFParts
{

class PairCorrelation : public SplinedFunction
{
public:
	PairCorrelation();
	PairCorrelation(int splineOrder_);

	int numOfParticlePairs;
	double numOfParticlePairs_Inv;
	int FirstParticleType;
	int SecondParticleType;

	bool useContactInteractionBoundaryCondition;
	bool useContactInteractionBoundaryConditionExp;
	double gamma;
	double contactLengthscale;

	double contactSum;
	double contactSumNew;
	vector<vector<double> > contactSumD;
	vector<double> contactSumD2;
	double contactOld;
	double contactNew;

public:
	void Init();
};

}
