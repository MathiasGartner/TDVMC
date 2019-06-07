#pragma once

#include <vector>

using namespace std;

class ParticleProperties
{
public:
	double mass;
	double hbarOver2m;

	double densityProfileMaxDistance;
	int numOfDensityProfileValues;
	vector<double> densityProfileBins;
	double densityProfileBinInterval;
	int densityProfileBin;
	double densityProfileNodePointSpacing;

public:
	ParticleProperties();
};
