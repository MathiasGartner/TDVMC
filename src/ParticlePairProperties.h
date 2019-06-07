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
