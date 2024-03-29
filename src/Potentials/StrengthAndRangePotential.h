/*-----------------------------------------------------------------------------
 *
 * 		Name:			StrengthAndRangePotential.h
 * 		Author:			Mathias Gartner
 * 		Description:	Prototype for an interaction potential characterized
 * 						by a strength parameter and a range parameter.
 *
 * --------------------------------------------------------------------------*/

#pragma once

#include "IPairPotential.h"

namespace Potentials
{

class StrengthAndRangePotential: public IPairPotential
{
public:
	double strength;
	double range;

public:
	StrengthAndRangePotential()
	{
		this->strength = 1.0;
		this->range = 1.0;
	}

	StrengthAndRangePotential(double strength_, double range_)
	{
		this->strength = strength_;
		this->range = range_;
	}
};

}
