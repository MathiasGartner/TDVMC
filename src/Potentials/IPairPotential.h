/*-----------------------------------------------------------------------------
 *
 * 		Name:			IPairPotential.h
 * 		Author:			Mathias Gartner
 * 		Description:	Interface for a distance dependent interaction
 * 						potential.
 *
 * --------------------------------------------------------------------------*/

#pragma once

#include <cmath>
#include <string>
#include <vector>

using namespace std;

namespace Potentials
{

class IPairPotential
{
public:
	string description;

public:
	IPairPotential()
	{
	}

	virtual ~IPairPotential()
	{
	}

	virtual double GetPotential(double r) = 0;
};

}

