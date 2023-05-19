/*-----------------------------------------------------------------------------
 *
 * 		Name:			ObservableVsOnGridWithScaling.h
 * 		Author:			Mathias Gartner
 * 		Description:	Vector valued observable on an underlying Grid
 * 						implementing volume-scaling in 1, 2 and 3 dimensions
 *
 * --------------------------------------------------------------------------*/

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
