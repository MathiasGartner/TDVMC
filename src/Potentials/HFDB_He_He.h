/*-----------------------------------------------------------------------------
 *
 * 		Name:			HFDB_He_He.h
 * 		Author:			Mathias Gartner
 * 		Description:	HFDB with parameters for Helium - Helium interaction
 * 						see
 * 							Aziz et al: Mol. Phys. 61, 1487 (1987)
 * 						see also HFDB.h
 *
 * --------------------------------------------------------------------------*/

#pragma once

#include "HFDB.h"

namespace Potentials
{

class HFDB_He_He: public HFDB
{
public:
	HFDB_He_He();
};

}
