#include "KTTY_He_Na.h"

using namespace std;

namespace Potentials
{

KTTY_He_Na::KTTY_He_Na()
{
	description = "KTTY_He_Na";

	d = 2.218564;
	b1 = 1.00872;
	b2 = 0.00399053;
	c6 = 23.768;
	c8 = 1307.6;
	c10 = 94563.2;

	c12 = pow(c10 / c8, 3.0) * c6;
	c14 = pow(c12 / c10, 3.0) * c8;
	c16 = pow(c14 / c12, 3.0) * c10;
}

}
