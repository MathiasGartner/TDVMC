/*-----------------------------------------------------------------------------
 *
 * 		Name:			ObservableV.h
 * 		Author:			Mathias Gartner
 * 		Description:	Vector valued observable
 *
 * --------------------------------------------------------------------------*/

#pragma once

#include "IObservable.h"

#include <vector>

using namespace std;

namespace Observables
{

class ObservableV : public IObservable
{
public:
	vector<double> values;

public:
	ObservableV();

	virtual ObservableV* Clone() const override;

	void ClearValues() override;
	void ApplySquareRoot() override;

	void Init(int size, string name);

	IObservable& operator+=(const IObservable& oc) override;
	IObservable& operator-=(const IObservable& oc) override;
	IObservable& operator*=(double d) override;
	IObservable& operator*=(const IObservable& oc) override;
	IObservable& operator/=(double d) override;
};

}
