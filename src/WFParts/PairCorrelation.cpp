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

	this->useContactInteractionBoundaryCondition = false;
	this->useContactInteractionBoundaryConditionExp = false;
	this->gamma = 0.0;
	this->contactLengthscale = 1.0;

	this->contactSum = 0.0;
	this->contactSumNew = 0.0;
	this->contactOld = 0;
	this->contactNew = 0;
}

void PairCorrelation::Init()
{
	SplinedFunction::Init();

	InitVector(this->contactSumD, N, DIM, 0.0);
	InitVector(this->contactSumD2, N, 0.0);
}

}
