#pragma once

#include "IPairPotential.h"

namespace Potentials
{

class LJ_He_He: public IPairPotential
{
public:
	double sigma;
	double eps;

public:
	LJ_He_He();

	double GetPotential(double r);
};

}
