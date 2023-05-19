/*-----------------------------------------------------------------------------
 *
 * 		Name:			ParticleProperties.h
 * 		Author:			Mathias Gartner
 * 		Description:	Holds properties for specific particle species.
 *
 * --------------------------------------------------------------------------*/

#pragma once

#include <vector>

using namespace std;

class ParticleProperties
{
public:
	double mass;
	double hbarOver2m;
	double kineticEnergyFactor;

	double densityProfileMaxDistance;
	int numOfDensityProfileValues;
	vector<double> densityProfileBins;
	double densityProfileBinInterval;
	int densityProfileBin;
	double densityProfileNodePointSpacing;

public:
	ParticleProperties();
};
