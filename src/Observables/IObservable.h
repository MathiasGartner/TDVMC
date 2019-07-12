#pragma once

#include <string>

using namespace std;

namespace Observables
{

class IObservable
{
public:
	string name;

public:
	IObservable()
	{
	}

	virtual ~IObservable()
	{
	}

	virtual IObservable* Clone() const = 0;

	virtual void ClearValues() = 0;

	virtual IObservable& operator+=(const IObservable& oc) = 0;
	virtual IObservable& operator-=(const IObservable& oc) = 0;
	virtual IObservable& operator*=(double d) = 0;
	virtual IObservable& operator/=(double d) = 0;
};

//template<typename Derived>
//class IObservable: public IIObservable
//{
//public:
//	virtual IIObservable *clone() const
//	{
//		return new Derived(static_cast<Derived const&>(*this));
//	}
//};

}

