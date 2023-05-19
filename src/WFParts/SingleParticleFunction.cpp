/*-----------------------------------------------------------------------------
 *
 * 		Name:			SingleParticleFunction.cpp
 * 		Author:			Mathias Gartner
 * 		Description:	see SingleParticleFunction.h
 *
 * --------------------------------------------------------------------------*/

#include "SingleParticleFunction.h"

namespace WFParts
{

SingleParticleFunction::SingleParticleFunction() : SplinedFunction(3)
{
}

SingleParticleFunction::SingleParticleFunction(int splineOrder_) : SplinedFunction(splineOrder_)
{
}

}
