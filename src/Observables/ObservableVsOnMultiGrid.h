/*-----------------------------------------------------------------------------
 *
 * 		Name:			ObservableVsOnMultiGrid.h
 * 		Author:			Mathias Gartner
 * 		Description:	Multi-dimensional observables on an underlying
 * 						Grid for each dimension. Used eg. for the pair density
 *
 * --------------------------------------------------------------------------*/

#pragma once

#include "Grid.h"
#include "IObservable.h"
#include "ObservableV.h"

#include <string>
#include <vector>

using namespace std;

namespace Observables
{

class ObservableVsOnMultiGrid : public IObservable
{
public:
	vector<ObservableV> observablesV;

	vector<Grid> grids;
	int totalGridPoints;

public:
	ObservableVsOnMultiGrid();

	virtual ObservableVsOnMultiGrid* Clone() const override;

	void ClearValues() override;
	void ApplySquareRoot() override;

	void InitGrid(vector<vector<double> > gridProperties);
	void InitObservables(vector<string> names);
	virtual void AddToHistogram(int observableIndex, vector<double> gridValues, double value);

	IObservable& operator+=(const IObservable& oc) override;
	IObservable& operator-=(const IObservable& oc) override;
	IObservable& operator*=(double d) override;
	IObservable& operator*=(const IObservable& oc) override;
	IObservable& operator/=(double d) override;

	vector<string> GridNames();
};

}
