#include "ObservableCollection.h"

namespace Observables
{

ObservableCollection::ObservableCollection()
{
}

ObservableCollection::~ObservableCollection()
{
}

ObservableCollection ObservableCollection::Clone() const
{
	ObservableCollection oc;
	for (auto& obs : this->observables)
	{
		oc.Add(obs->Clone());
	}
	return oc;
}

void ObservableCollection::Destroy()
{
	for (auto& obs : this->observables)
	{
		delete obs;
	}
}

void ObservableCollection::Add(IObservable* obs)
{
	this->observables.push_back(obs);
}

void ObservableCollection::ClearValues()
{
	for (auto& obs : this->observables)
	{
		obs->ClearValues();
	}
}

void ObservableCollection::ApplySquareRoot()
{
	for (auto& obs : this->observables)
	{
		obs->ApplySquareRoot();
	}
}

ObservableCollection& ObservableCollection::operator+=(const ObservableCollection& oc)
{
	for (unsigned int i = 0; i < this->observables.size(); i++)
	{
		*(this->observables[i]) += *(oc.observables[i]);
	}
	return *this;
}

ObservableCollection& ObservableCollection::operator-=(const ObservableCollection& oc)
{
	for (unsigned int i = 0; i < this->observables.size(); i++)
	{
		*(this->observables[i]) -= *(oc.observables[i]);
	}
	return *this;
}

ObservableCollection& ObservableCollection::operator*=(double d)
{
	for (unsigned int i = 0; i < this->observables.size(); i++)
	{
		*(this->observables[i]) *= d;
	}
	return *this;
}

ObservableCollection& ObservableCollection::operator*=(ObservableCollection& oc)
{
	for (unsigned int i = 0; i < this->observables.size(); i++)
	{
		*(this->observables[i]) *= *(oc.observables[i]);
	}
	return *this;
}

ObservableCollection& ObservableCollection::operator/=(double d)
{
	for (unsigned int i = 0; i < this->observables.size(); i++)
	{
		*(this->observables[i]) /= d;
	}
	return *this;
}

}
