/*-----------------------------------------------------------------------------
 *
 * 		Name:			LJ.h
 * 		Author:			Mathias Gartner
 * 		Description:	Lennard-Jones potential.
 *
 * --------------------------------------------------------------------------*/

#pragma once

#include "IPairPotential.h"

namespace Potentials
{

class LJ: public IPairPotential
{
public:
	double sigma;
	double eps;

public:
	LJ();

	double GetPotential(double r) override;
};

}
