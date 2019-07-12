#pragma once

#include "ObservableVsOnGrid.h"

using namespace std;

namespace Observables
{

class ObservableVsOnGridWithScaling : public ObservableVsOnGrid
{
public:
	vector<double> scalingGrid;

public:
	ObservableVsOnGridWithScaling();

	virtual ObservableVsOnGridWithScaling* Clone() const override;

	void InitScaling();
	void AddToHistogram(int observableIndex, double gridValue, double value) override;
};

}
