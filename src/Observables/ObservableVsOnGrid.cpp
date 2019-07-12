#include "ObservableVsOnGrid.h"

#include "../Utils.h"

namespace Observables
{

ObservableVsOnGrid::ObservableVsOnGrid()
{
	gridName = "grid";
	gridMin = 0.0;
	gridMax = 0.0;
	gridCount = 0;
	gridSpacing = 0.0;
}

ObservableVsOnGrid* ObservableVsOnGrid::Clone() const
{
	return new ObservableVsOnGrid(*this);
}

void ObservableVsOnGrid::ClearValues()
{
	for (auto& v : this->observablesV)
	{
		v.ClearValues();
	}
}

void ObservableVsOnGrid::InitGrid(double min, double max, double spacing)
{
	gridMin = min;
	gridMax = max;
	gridSpacing = spacing;
	gridCount = (gridMax - gridMin) / gridSpacing;
	InitVector(grid, gridCount, 0.0);

	double val = gridMin;
	for (unsigned int i = 0; i < grid.size(); i++)
	{
		grid[i] = val;
		val += gridSpacing;

	}
}

void ObservableVsOnGrid::InitObservables(vector<string> names)
{
	this->observablesV.resize(names.size());
	for(unsigned int i = 0; i < names.size(); i++)
	{
		this->observablesV[i].Init(this->grid.size(), names[i]);
	}
}

void ObservableVsOnGrid::AddToHistogram(int observableIndex, double gridValue, double value)
{
	//INFO: only for grid that starts at zero. otherwise: (gridValue - this->gridMin) has to be used.
	double interval;
	int bin;
	interval = gridValue / this->gridSpacing;
	bin = floor(interval);
	this->observablesV[observableIndex].values[bin] += value;
}

IObservable& ObservableVsOnGrid::operator+=(const IObservable& oc)
{
	for (unsigned int i = 0; i < this->observablesV.size(); i++)
	{
		this->observablesV[i] += dynamic_cast<const ObservableVsOnGrid&>(oc).observablesV[i];
	}
	return *this;
}

IObservable& ObservableVsOnGrid::operator-=(const IObservable& oc)
{
	for (unsigned int i = 0; i < this->observablesV.size(); i++)
	{
		this->observablesV[i] -= dynamic_cast<const ObservableVsOnGrid&>(oc).observablesV[i];
	}
	return *this;
}

IObservable& ObservableVsOnGrid::operator*=(double d)
{
	for (unsigned int i = 0; i < this->observablesV.size(); i++)
	{
		this->observablesV[i] *= d;
	}
	return *this;
}

IObservable& ObservableVsOnGrid::operator/=(double d)
{
	for (unsigned int i = 0; i < this->observablesV.size(); i++)
	{
		this->observablesV[i] /= d;
	}
	return *this;
}

}
