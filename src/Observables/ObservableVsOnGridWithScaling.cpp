#include "ObservableVsOnGridWithScaling.h"

#include "../Utils.h"

namespace Observables
{

ObservableVsOnGridWithScaling::ObservableVsOnGridWithScaling()
{

}

ObservableVsOnGridWithScaling* ObservableVsOnGridWithScaling::Clone() const
{
	return new ObservableVsOnGridWithScaling(*this);
}

//INFO: only for r^2 volume observables in 3D...
void ObservableVsOnGridWithScaling::InitScaling()
{
	InitVector(scalingGrid, gridCount, 0.0);
	for (int i = 0; i < gridCount; i++)
	{
		scalingGrid[i] = 4.0 * M_PI * pow(gridSpacing * (i + 1), 3.0) / 3.0;
	}
	for (int i = gridCount - 1; i > 0; i--)
	{
		scalingGrid[i] = scalingGrid[i] - scalingGrid[i - 1];
	}
}

//TODO: properly override using virtual function in base class...
void ObservableVsOnGridWithScaling::AddToHistogram(int observableIndex, double gridValue, double value)
{
	//INFO: only for grid that starts at zero. otherwise: (gridValue - this->gridMin) has to be used.
	double interval;
	int bin;
	interval = gridValue / this->gridSpacing;
	bin = floor(interval);
	this->observablesV[observableIndex].values[bin] += value / this->scalingGrid[bin];
}

}
