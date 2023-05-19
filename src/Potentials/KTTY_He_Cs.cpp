/*-----------------------------------------------------------------------------
 *
 * 		Name:			KTTY_He_Cs.cpp
 * 		Author:			Mathias Gartner
 * 		Description:	see KTTY_He_Cs.h
 *
 * --------------------------------------------------------------------------*/

#include "KTTY_He_Cs.h"

using namespace std;

namespace Potentials
{

KTTY_He_Cs::KTTY_He_Cs()
{
	description = "KTTY_He_Cs";

	d = 1.224951;
	b1 = 0.782095;
	b2 = 0.00513175;
	c6 = 41.417;
	c8 = 3903.4;
	c10 = 453443.0;

	c12 = pow(c10 / c8, 3.0) * c6;
	c14 = pow(c12 / c10, 3.0) * c8;
	c16 = pow(c14 / c12, 3.0) * c10;
}

}
