/*-----------------------------------------------------------------------------
 *
 * 		Name:			SquareWell.h
 * 		Author:			Mathias Gartner
 * 		Description:	Rectangular shaped potential.
 *
 * --------------------------------------------------------------------------*/

#pragma once

#include "StrengthAndRangePotential.h"

namespace Potentials
{

class SquareWell: public StrengthAndRangePotential
{

public:
	SquareWell();
	SquareWell(double strength_, double range_);

	double GetPotential(double r) override;
};

}
