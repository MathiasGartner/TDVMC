#pragma once

#include "IPairPotential.h"

namespace Potentials
{

class None: public IPairPotential
{
public:

public:
	None();

	double GetPotential(double r) override;
};

}
