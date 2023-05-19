/*-----------------------------------------------------------------------------
 *
 * 		Name:			ParticlePairProperties.h
 * 		Author:			Mathias Gartner
 * 		Description:	Holds properties for a pair of particles of different
 * 						(or same) particle type.
 *
 * --------------------------------------------------------------------------*/

#pragma once

#include "Potentials/IPairPotential.h"

#include <vector>

using namespace std;

class ParticlePairProperties
{
public:
	Potentials::IPairPotential* potential;

public:
	ParticlePairProperties();
};
