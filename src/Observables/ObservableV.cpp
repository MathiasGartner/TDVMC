/*-----------------------------------------------------------------------------
 *
 * 		Name:			ObservableV.cpp
 * 		Author:			Mathias Gartner
 * 		Description:	see ObservableV.h
 *
 * --------------------------------------------------------------------------*/

#include "ObservableV.h"

#include "../Utils.h"

namespace Observables
{

ObservableV::ObservableV()
{

}

ObservableV* ObservableV::Clone() const
{
	return new ObservableV(*this);
}

void ObservableV::ClearValues()
{
	ClearVector(this->values);
}

void ObservableV::ApplySquareRoot()
{
	for (int i = 0; i < this->values.size(); i++)
	{
		this->values[i] = sqrt(this->values[i]);
	}
}

void ObservableV::Init(int size, string name)
{
	InitVector(this->values, size, 0.0);
	this->name = name;
}

IObservable& ObservableV::operator+=(const IObservable& oc)
{
	this->values += dynamic_cast<const ObservableV&>(oc).values;
	return *this;
}

IObservable& ObservableV::operator-=(const IObservable& oc)
{
	this->values -= dynamic_cast<const ObservableV&>(oc).values;
	return *this;
}

IObservable& ObservableV::operator*=(double d)
{
	this->values *= d;
	return *this;
}

IObservable& ObservableV::operator*=(const IObservable& oc)
{
	this->values *= dynamic_cast<const ObservableV&>(oc).values;
	return *this;
}

IObservable& ObservableV::operator/=(double d)
{
	this->values /= d;
	return *this;
}

}
