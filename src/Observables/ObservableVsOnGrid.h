#pragma once

#include "ObservableV.h"
#include "IObservable.h"

#include <string>
#include <vector>

using namespace std;

namespace Observables
{

class ObservableVsOnGrid : public IObservable
{
public:
	vector<ObservableV> observablesV;

	string gridName;
	double gridMin;
	double gridMax;
	double gridCount;
	double gridSpacing;
	vector<double> grid;

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
