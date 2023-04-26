#include "ObservableVsOnGrid.h"

#include "../Utils.h"

namespace Observables
{

ObservableVsOnGrid::ObservableVsOnGrid()
{
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

void ObservableVsOnGrid::ApplySquareRoot()
{
	for (auto& v : this->observablesV)
	{
		v.ApplySquareRoot();
	}
}

void ObservableVsOnGrid::InitGrid(double min, double max, double spacing)
{
	this->grid.Init(min, max, spacing);
}

void ObservableVsOnGrid::InitGrid(vector<double> g)
{
	this->grid.Init(g);
}

void ObservableVsOnGrid::InitObservables(vector<string> names)
{
	this->observablesV.resize(names.size());
	for(unsigned int i = 0; i < names.size(); i++)
	{
		this->observablesV[i].Init(this->grid.count, names[i]);
	}
}

void ObservableVsOnGrid::AddToHistogram(int observableIndex, double gridValue, double value)
{
	int bin;
	bin = this->grid.GetIndex(gridValue);
	this->observablesV[observableIndex].values[bin] += value;
}

void ObservableVsOnGrid::SetValueAtGridIndex(int observableIndex, double gridIndex, double value)
{
	this->observablesV[observableIndex].values[gridIndex] = value;
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

IObservable& ObservableVsOnGrid::operator*=(const IObservable& oc)
{
	for (unsigned int i = 0; i < this->observablesV.size(); i++)
	{
		this->observablesV[i] *= dynamic_cast<const ObservableVsOnGrid&>(oc).observablesV[i];
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
