/*-----------------------------------------------------------------------------
 *
 * 		Name:			LJ_XSplit.cpp
 * 		Author:			Mathias Gartner
 * 		Description:	see LJ_XSplit.h
 *
 * --------------------------------------------------------------------------*/

#include "LJ_XSplit.h"

using namespace std;

namespace Potentials
{

LJ_XSplit::LJ_XSplit()
{
	description = "LJ_XSplit";

	//X_Cluster Split
	sigma = 4.0;
	eps = 3.56;
}

}
