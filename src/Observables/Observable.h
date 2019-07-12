#pragma once

#include "IObservable.h"

using namespace std;

namespace Observables
{

class Observable : public IObservable
{
public:
	double value;

public:
	Observable();

	virtual Observable* Clone() const override;

	void ClearValues() override;

	IObservable& operator+=(const IObservable& oc) override;
	IObservable& operator-=(const IObservable& oc) override;
	IObservable& operator*=(double d) override;
	IObservable& operator/=(double d) override;
};

}
