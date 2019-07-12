#include "LJ.h"

using namespace std;

namespace Potentials
{

LJ::LJ()
{
	description = "LJ";

	sigma = 1.0;
	eps = 1.0;
}

double LJ::GetPotential(double r)
{
	double LJ;
	double sigma_r_6 = pow(sigma / r, 6);
	LJ = 4.0 * eps * sigma_r_6 * (sigma_r_6 - 1.0);
	return LJ;
}

}
