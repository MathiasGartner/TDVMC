/*-----------------------------------------------------------------------------
 *
 * 		Name:			LJ_XSplit.h
 * 		Author:			Mathias Gartner
 * 		Description:	Lennard-Jones potential with parameter values used
 * 						in test systems by Petar Stipanovic in Split, Croatia
 *
 * --------------------------------------------------------------------------*/

#pragma once

#include "LJ.h"

namespace Potentials
{

class LJ_XSplit: public LJ
{
public:
	LJ_XSplit();
};

}
