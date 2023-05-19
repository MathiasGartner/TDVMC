/*-----------------------------------------------------------------------------
 *
 * 		Name:			ObservableCollection.h
 * 		Author:			Mathias Gartner
 * 		Description:	Holds a list of elements implementing the IObservable
 * 						interface. Every IPhysicalSystem should hold an
 * 						ObservableCollection to store the data.
 *
 * --------------------------------------------------------------------------*/

#pragma once

#include "IObservable.h"

#include <vector>

using namespace std;

namespace Observables
{

class ObservableCollection
{
public:
	//TODO: should be private...
	vector<IObservable*> observables;

public:
	ObservableCollection();
	~ObservableCollection();

	ObservableCollection Clone() const;

	void Destroy();

	void Add(IObservable* obs);

	void ClearValues();
	void ApplySquareRoot();

	ObservableCollection& operator+=(const ObservableCollection& oc);
	ObservableCollection& operator-=(const ObservableCollection& oc);
	ObservableCollection& operator*=(double d);
	ObservableCollection& operator*=(ObservableCollection& oc);
	ObservableCollection& operator/=(double d);

};

}
