#pragma once

#include "Grid.h"
#include "IObservable.h"
#include "ObservableV.h"

#include <string>
#include <vector>

using namespace std;

namespace Observables
{

class ObservableVsOnGrid : public IObservable
{
public:
	vector<ObservableV> observablesV;

	Grid grid;

public:
	ObservableVsOnGrid();

	virtual ObservableVsOnGrid* Clone() const override;

	void ClearValues() override;

	void InitGrid(double min, double max, double spacing);
	void InitObservables(vector<string> names);
	virtual void AddToHistogram(int observableIndex, double gridValue, double value);

	IObservable& operator+=(const IObservable& oc) override;
	IObservable& operator-=(const IObservable& oc) override;
	IObservable& operator*=(double d) override;
	IObservable& operator/=(double d) override;
};

}