#include "PairCorrelation.h"

namespace WFParts
{

PairCorrelation::PairCorrelation() : SplinedFunction(3)
{
	this->FirstParticleType = -1;
	this->SecondParticleType = -1;
}

PairCorrelation::PairCorrelation(int splineOrder_) : SplinedFunction(splineOrder_)
{
	this->FirstParticleType = -1;
	this->SecondParticleType = -1;
}

}
