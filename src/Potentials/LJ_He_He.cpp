#include "LJ_He_He.h"

using namespace std;

namespace Potentials
{

LJ_He_He::LJ_He_He()
{
	sigma = 4.0;
	eps = 3.56;
}

double LJ_He_He::GetPotential(double r)
{
	double LJ;
	double sigma_r_6 = pow(sigma / r, 6);
	LJ = 4.0 * eps * sigma_r_6 * (sigma_r_6 - 1.0);
	return LJ;
}

}
