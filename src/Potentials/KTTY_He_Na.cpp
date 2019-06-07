#include "KTTY_He_Na.h"

using namespace std;

namespace Potentials
{

KTTY_He_Na::KTTY_He_Na()
{
	epsil = 3.1577504e8;
	d = 2.218564;
	b1 = 1.00872;
	b2 = 0.00399053;
	c6 = 23.768;
	c8 = 1307.6;
	c10 = 94563.2;

	c12 = pow(c10 / c8, 3.0) * c6;
	c14 = pow(c12 / c10, 3.0) * c8;
	c16 = pow(c14 / c12, 3.0) * c10;

	fak2 = 2.0;
	fak3 = fak2 * 3.0;
	fak4 = fak3 * 4.0;
	fak5 = fak4 * 5.0;
	fak6 = fak5 * 6.0;
	fak7 = fak6 * 7.0;
	fak8 = fak7 * 8.0;
	fak9 = fak8 * 9.0;
	fak10 = fak9 * 10.0;
	fak11 = fak10 * 11.0;
	fak12 = fak11 * 12.0;
	fak13 = fak12 * 13.0;
	fak14 = fak13 * 14.0;
	fak15 = fak14 * 15.0;
	fak16 = fak15 * 16.0;
}

double KTTY_He_Na::GetPotential(double r)
{

	double x = r / 0.52917721092;
	double x2 = x * x;
	double xminus2 = 1.0 / x2;
	double xminus6 = xminus2 * xminus2 * xminus2;
	double xminus8 = xminus6 * xminus2;
	double xminus10 = xminus8 * xminus2;
	double xminus12 = xminus10 * xminus2;
	double xminus14 = xminus12 * xminus2;
	double xminus16 = xminus14 * xminus2;

	double bet = b1 * x + b2 * x * x;

	double vrep = d * exp(-bet);
	double br = (b1 + 2.0 * b2 * x) * x;
	double exbr = exp(-br);
	double br2 = br * br;
	double br3 = br2 * br;
	double br4 = br3 * br;
	double br5 = br4 * br;
	double br6 = br5 * br;
	double br7 = br6 * br;
	double br8 = br7 * br;
	double br9 = br8 * br;
	double br10 = br9 * br;
	double br11 = br10 * br;
	double br12 = br11 * br;
	double br13 = br12 * br;
	double br14 = br13 * br;
	double br15 = br14 * br;
	double br16 = br15 * br;
	double f6 = 1.0 - exbr * (1.0 + br + br2 / fak2 + br3 / fak3 + br4 / fak4 + br5 / fak5 + br6 / fak6);
	double f8 = f6 - exbr * (br7 / fak7 + br8 / fak8);
	double f10 = f8 - exbr * (br9 / fak9 + br10 / fak10);
	double f12 = f10 - exbr * (br11 / fak11 + br12 / fak12);
	double f14 = f12 - exbr * (br13 / fak13 + br14 / fak14);
	double f16 = f14 - exbr * (br15 / fak15 + br16 / fak16);

	double vatr = f6 * c6 * xminus6 + f8 * c8 * xminus8 + f10 * c10 * xminus10 + f12 * c12 * xminus12 + f14 * c14 * xminus14 + f16 * c16 * xminus16;
	double fpot = vrep - vatr;
	double fpotHeA = epsil * fpot;

	return fpotHeA; //TODO: check if in mK or K!
}

}
