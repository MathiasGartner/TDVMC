/*-----------------------------------------------------------------------------
 *
 * 		Name:			Gauss.h
 * 		Author:			Mathias Gartner
 * 		Description:	Gaussian shaped potential.
 *
 * --------------------------------------------------------------------------*/

#pragma once

#include "StrengthAndRangePotential.h"

namespace Potentials
{

class Gauss: public StrengthAndRangePotential
{

public:
	Gauss();
	Gauss(double strength_, double range_);

	double GetPotential(double r) override;
};

}
