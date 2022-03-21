#pragma once

#include "SplinedFunction.h"

#include <vector>

using namespace std;

namespace WFParts
{

class PairCorrelation : public SplinedFunction
{
public:
	PairCorrelation();
	PairCorrelation(int splineOrder_);

	int numOfParticlePairs;
	int numOfParticlePairs_Inv;
	int FirstParticleType;
	int SecondParticleType;
};

}
