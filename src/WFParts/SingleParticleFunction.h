#pragma once

#include "SplinedFunction.h"

#include <vector>

using namespace std;

namespace WFParts
{

class SingleParticleFunction : public SplinedFunction
{
public:
	SingleParticleFunction();
	SingleParticleFunction(int splineOrder_);
};

}
