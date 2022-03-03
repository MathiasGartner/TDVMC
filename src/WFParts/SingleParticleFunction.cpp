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
