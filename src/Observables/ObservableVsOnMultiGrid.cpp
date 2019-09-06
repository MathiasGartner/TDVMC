#include "ObservableVsOnMultiGrid.h"

#include "../Utils.h"

namespace Observables
{

ObservableVsOnMultiGrid::ObservableVsOnMultiGrid()
{
	totalGridPoints = 0;
}

ObservableVsOnMultiGrid* ObservableVsOnMultiGrid::Clone() const
{
	return new ObservableVsOnMultiGrid(*this);
}

void ObservableVsOnMultiGrid::ClearValues()
{
	for (auto& v : this->observablesV)
	{
		v.ClearValues();
	}
}

void ObservableVsOnMultiGrid::InitGrid(vector<vector<double> > gridProperties)
{
	totalGridPoints = 1;
	this->grids.resize(gridProperties.size());
	int i = 0;
	for (auto& g : grids)
	{
		g.min = gridProperties[i][0];
		g.max = gridProperties[i][1];
		g.spacing = gridProperties[i][2];
		g.count = (g.max - g.min) / g.spacing;
		InitVector(g.grid, g.count, 0.0);

		double val = g.min;
		for (int j = 0; j < g.count; j++)
		{
			g.grid[j] = val;
			val += g.spacing;
		}
		totalGridPoints *= g.count;
		i++;
	}
}

void ObservableVsOnMultiGrid::InitObservables(vector<string> names)
{
	this->observablesV.resize(names.size());
	for(unsigned int i = 0; i < names.size(); i++)
	{
		this->observablesV[i].Init(this->totalGridPoints, names[i]);
	}
}

void ObservableVsOnMultiGrid::AddToHistogram(int observableIndex, vector<double> gridValues, double value)
{
	int indexForObservableValue = 0;
	int offset = 1;
	for (unsigned int i = 0; i < this->grids.size(); i++)
	{
		indexForObservableValue += this->grids[i].GetIndex(gridValues[i]) * offset;
		offset *= this->grids[i].count;
	}
	this->observablesV[observableIndex].values[indexForObservableValue] += value;
}

IObservable& ObservableVsOnMultiGrid::operator+=(const IObservable& oc)
{
	for (unsigned int i = 0; i < this->observablesV.size(); i++)
	{
		this->observablesV[i] += dynamic_cast<const ObservableVsOnMultiGrid&>(oc).observablesV[i];
	}
	return *this;
}

IObservable& ObservableVsOnMultiGrid::operator-=(const IObservable& oc)
{
	for (unsigned int i = 0; i < this->observablesV.size(); i++)
	{
		this->observablesV[i] -= dynamic_cast<const ObservableVsOnMultiGrid&>(oc).observablesV[i];
	}
	return *this;
}

IObservable& ObservableVsOnMultiGrid::operator*=(double d)
{
	for (unsigned int i = 0; i < this->observablesV.size(); i++)
	{
		this->observablesV[i] *= d;
	}
	return *this;
}

IObservable& ObservableVsOnMultiGrid::operator/=(double d)
{
	for (unsigned int i = 0; i < this->observablesV.size(); i++)
	{
		this->observablesV[i] /= d;
	}
	return *this;
}

}
