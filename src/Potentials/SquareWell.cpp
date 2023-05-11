/*-----------------------------------------------------------------------------
 *
 * 		Name:			SquareWell.cpp
 * 		Author:			Mathias Gartner
 * 		Description:	see SquareWell.h
 *
 * --------------------------------------------------------------------------*/

#include "SquareWell.h"

using namespace std;

namespace Potentials
{

SquareWell::SquareWell() : StrengthAndRangePotential()
{
	this->description = "SquareWell";
}

SquareWell::SquareWell(double strength_, double range_) : StrengthAndRangePotential(strength_, range_)
{
	this->description = "SquareWell";
}

double SquareWell::GetPotential(double r)
{
	double pot = 0;
	if (r < this->range)
	{
		pot = this->strength;
	}
	return pot;
}

}
