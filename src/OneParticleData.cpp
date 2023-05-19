/*-----------------------------------------------------------------------------
 *
 * 		Name:			OneParticleData.cpp
 * 		Author:			Mathias Gartner
 * 		Description:	see OneParticleData.h
 *
 * --------------------------------------------------------------------------*/

#include "OneParticleData.h"

OneParticleData::OneParticleData()
{
	kineticSumR1 = 0;
	kineticSumI1 = 0;
	kineticSumR1I1 = 0;
	kineticSumR2 = 0;
	kineticSumI2 = 0;

	vecKineticSumR1.resize(DIM);
	vecKineticSumI1.resize(DIM);
}
