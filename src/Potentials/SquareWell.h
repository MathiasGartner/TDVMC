#pragma once

#include "IPairPotential.h"

namespace Potentials
{

class SquareWell: public IPairPotential
{
public:
	double strength;
	double width;

public:
	SquareWell();

	double GetPotential(double r) override;
};

}
