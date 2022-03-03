#include "SquareWell.h"

using namespace std;

namespace Potentials
{

SquareWell::SquareWell()
{
	this->description = "SquareWell";

	this->strength = 1.0;
	this->width = 1.0;
}

double SquareWell::GetPotential(double r)
{
	double pot = 0;
	if (r < width)
	{
		pot = strength;
	}
	return pot;
}

}
