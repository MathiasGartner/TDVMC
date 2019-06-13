#pragma once

#include "../Utils.h"

#include "HFDB_He_He.h"
#include "KTTY_He_Na.h"
#include "LJ_He_He.h"

#include <string>
#include <vector>

using namespace std;

namespace Potentials
{

void SavePotentialValues(IPairPotential& pot, double rmin, double rmax, double step, string filename)
{
	vector<vector<double> > values;
	for (double r = rmin; r <= rmax; r += step)
	{
		values.push_back( { r, pot.GetPotential(r) });
	}
	WriteDataToFile(values, filename, pot.description);
}

void SaveAllPotentialValues(string path)
{
	HFDB_He_He hfdb_he_he;
	SavePotentialValues(hfdb_he_he, 0.1, 20.0, 0.1, path + "HFDB_He_He");

	KTTY_He_Na ktty_he_na;
	SavePotentialValues(ktty_he_na, 0.1, 20.0, 0.1, path + "KTTY_He_Na");

	LJ_He_He lj_he_he;
	SavePotentialValues(lj_he_he, 0.1, 20.0, 0.1, path + "LJ_He_He");
}

}
