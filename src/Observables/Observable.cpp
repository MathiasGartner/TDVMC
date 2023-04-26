#include "Observable.h"

namespace Observables
{

Observable::Observable()
{
	value = 0.0;
}

Observable* Observable::Clone() const
{
	return new Observable(*this);
}

void Observable::ClearValues()
{
	this->value = 0.0;
}

void Observable::ApplySquareRoot()
{
	this->value = sqrt(this->value);
}

IObservable& Observable::operator+=(const IObservable& oc)
{
	this->value += dynamic_cast<const Observable&>(oc).value;
	return *this;
}

IObservable& Observable::operator-=(const IObservable& oc)
{
	this->value -= dynamic_cast<const Observable&>(oc).value;
	return *this;
}

IObservable& Observable::operator*=(double d)
{
	this->value *= d;
	return *this;
}

IObservable& Observable::operator*=(const IObservable& oc)
{
	this->value *= dynamic_cast<const Observable&>(oc).value;
	return *this;
}

IObservable& Observable::operator/=(double d)
{
	this->value /= d;
	return *this;
}

}
