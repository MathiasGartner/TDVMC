/*-----------------------------------------------------------------------------
 *
 * 		Name:			LJ_He_He.cpp
 * 		Author:			Mathias Gartner
 * 		Description:	see LJ_He_He.h
 *
 * --------------------------------------------------------------------------*/

#include "LJ_He_He.h"

using namespace std;

namespace Potentials
{

LJ_He_He::LJ_He_He()
{
	description = "LJ_He_He";

	//He4
	sigma = 2.556;
	eps = 10.22;
}
}
