#pragma once

#include "../Utils.h"

#include "HFDB_He_He.h"
#include "KTTY_He_Cs.h"
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
	for (double r = rmin; r <= rmax + step / 2.0; r += step)
	{
		values.push_back( { r, pot.GetPotential(r) });
	}
	WriteDataToFile(values, filename, { "r", pot.description });
}

void SaveAllPotentialValues(string path)
{
	double min = 1.0;
	double max = 20.0;
	double step = 0.01;

	HFDB_He_He hfdb_he_he;
	SavePotentialValues(hfdb_he_he, min, max, step, path + "HFDB_He_He");

	KTTY_He_Cs ktty_he_cs;
	SavePotentialValues(ktty_he_cs, min, max, step, path + "KTTY_He_Cs");

	KTTY_He_Na ktty_he_na;
	SavePotentialValues(ktty_he_na, min, max, step, path + "KTTY_He_Na");

	LJ_He_He lj_he_he;
	SavePotentialValues(lj_he_he, min, max, step, path + "LJ_He_He");
}

}
