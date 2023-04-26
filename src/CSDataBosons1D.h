#pragma once

#include "ICorrelatedSamplingData.h"

#include "Observables/ObservableCollection.h"

#include <vector>

using namespace std;

class CSDataBosons1D: public ICorrelatedSamplingData
{
public:
	Observables::ObservableCollection data;
};
