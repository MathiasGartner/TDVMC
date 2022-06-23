#include "ParticleProperties.h"

ParticleProperties::ParticleProperties()
{
	mass = 0;
	hbarOver2m = 0;
	kineticEnergyFactor = 0.0;

	densityProfileMaxDistance = 0;
	numOfDensityProfileValues = 0;
	densityProfileBinInterval = 0;
	densityProfileBin = 0;
	densityProfileNodePointSpacing = 0;
}
