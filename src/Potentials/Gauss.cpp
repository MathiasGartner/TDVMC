#include "Gauss.h"

using namespace std;

namespace Potentials
{

Gauss::Gauss() : StrengthAndRangePotential(100.0, 0.1)
{
	this->description = "Gauss";
}

Gauss::Gauss(double strength_, double range_) : StrengthAndRangePotential(strength_, range_)
{
	this->description = "Gauss";
}

double Gauss::GetPotential(double r)
{
	double pot = 0;
	double r_a = r / this->range;
	pot = this->strength * exp(-(r_a * r_a) / 2.0);
	return pot;
}

}
