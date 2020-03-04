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
};

}
