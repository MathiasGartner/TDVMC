/*-----------------------------------------------------------------------------
 *
 * 		Name:			None.cpp
 * 		Author:			Mathias Gartner
 * 		Description:	see None.h
 *
 * --------------------------------------------------------------------------*/

#include "None.h"

using namespace std;

namespace Potentials
{

None::None()
{
	description = "None";
}

double None::GetPotential(double r)
{
	return 0.0;
}

}
