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
	if (DIM == 1)
	{
		InitVector(scalingGrid, grid.count, grid.spacing);
	}
	else
	{
		InitVector(scalingGrid, grid.count, 0.0);
		for (int i = 0; i < grid.count; i++)
		{
			if (DIM == 3)
			{
				scalingGrid[i] = 4.0 * M_PI * pow(grid.spacing * (i + 1), 3.0) / 3.0; //INFO: 3D sphere volume
			}
			else if (DIM == 2)
			{
				scalingGrid[i] = M_PI * pow(grid.spacing * (i + 1), 2.0); //INFO: 2D "sphere" volume
			}
		}
		for (int i = grid.count - 1; i > 0; i--)
		{
			scalingGrid[i] = scalingGrid[i] - scalingGrid[i - 1];
		}
	}
}

//TODO: properly override using virtual function in base class...
void ObservableVsOnGridWithScaling::AddToHistogram(int observableIndex, double gridValue, double value)
{
	int bin;
	bin = this->grid.GetIndex(gridValue);
	this->observablesV[observableIndex].values[bin] += value / this->scalingGrid[bin];
}

}
