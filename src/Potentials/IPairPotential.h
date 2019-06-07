#pragma once

#include <cmath>
#include <vector>

using namespace std;

namespace Potentials
{

class IPairPotential
{
public:

public:
	IPairPotential()
	{
	}

	virtual double GetPotential(double r) = 0;

	virtual ~IPairPotential()
	{
	}
};

}

