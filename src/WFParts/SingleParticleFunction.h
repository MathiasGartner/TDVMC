/*-----------------------------------------------------------------------------
 *
 * 		Name:			SplinedFunction.h
 * 		Author:			Mathias Gartner
 * 		Description:	Single particle function of the variational wavefunction
 *
 * --------------------------------------------------------------------------*/

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
