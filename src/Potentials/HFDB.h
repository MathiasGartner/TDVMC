/*-----------------------------------------------------------------------------
 *
 * 		Name:			HFDB.h
 * 		Author:			Mathias Gartner
 * 		Description:	HFDB potential, see:
 * 							Aziz et al: Mol. Phys. 61, 1487 (1987)
 *
 * --------------------------------------------------------------------------*/

#pragma once

#include "IPairPotential.h"

namespace Potentials
{

class HFDB: public IPairPotential
{
public:
	double epsil;
	double rm;
	double av;
	double alf;
	double bet;
	double dv;
	double c6;
	double c8;
	double c10;

public:
	HFDB();

	double GetPotential(double r) override;
};

}
