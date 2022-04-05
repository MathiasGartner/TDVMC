#include "PairCorrelation.h"

namespace WFParts
{

PairCorrelation::PairCorrelation() : PairCorrelation(3)
{
}

PairCorrelation::PairCorrelation(int splineOrder_) : SplinedFunction(splineOrder_)
{
	this->FirstParticleType = -1;
	this->SecondParticleType = -1;
	this->numOfParticlePairs = 0;
	this->numOfParticlePairs_Inv = 0.0;
}

}
