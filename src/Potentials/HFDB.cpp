#include "HFDB.h"

using namespace std;

namespace Potentials
{

HFDB::HFDB()
{
	description = "HFDB";

	epsil = 10.948;
	rm = 2.9630;
	av = 184431.01;
	alf = 10.43329537;
	bet = -2.27965105;
	dv = 1.4826;
	c6 = 1.36745214;
	c8 = 0.42123807;
	c10 = 0.17473318;
}

double HFDB::GetPotential(double r)
{
	double fpot = 0;

	double x = r / rm;
	double x2 = x * x;
	double xminus2 = 1.0 / x2;
	double xminus6 = xminus2 * xminus2 * xminus2;
	double xminus8 = xminus6 * xminus2;
	double xminus10 = xminus8 * xminus2;
	double f3 = c6 * xminus6 + c8 * xminus8 + c10 * xminus10;
	double f4 = av * exp(-alf * x + bet * x2);
	if (x >= dv)
	{
		fpot = f4 - f3;
	}
	else
	{
		double tmp = dv / x - 1.0;
		double f2 = exp(-(tmp * tmp));
		fpot = f4 - f3 * f2;
	}
	fpot = epsil * fpot;
	return fpot;
}

}
