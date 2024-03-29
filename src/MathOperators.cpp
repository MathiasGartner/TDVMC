/*-----------------------------------------------------------------------------
 *
 * 		Name:			MathOperators.h
 * 		Author:			Mathias Gartner
 * 		Description:	see MathOperators.h
 *
 * --------------------------------------------------------------------------*/

#include "MathOperators.h"

using namespace std;

vector<double> operator*(vector<int> v, double x)
{
	vector<double> result(v.size());
	for (unsigned int i = 0; i < v.size(); i++)
	{
		result[i] = v[i] * x;
	}
	return result;
}

vector<double>& operator+=(vector<double>& v1, const vector<double>& v2)
{
	for (unsigned int i = 0; i < v1.size(); i++)
	{
		v1[i] += v2[i];
	}
	return v1;
}

vector<double>& operator+=(vector<double>& v1, double d)
{
	for (unsigned int i = 0; i < v1.size(); i++)
	{
		v1[i] += d;
	}
	return v1;
}

vector<double>& operator-=(vector<double>& v1, const vector<double>& v2)
{
	for (unsigned int i = 0; i < v1.size(); i++)
	{
		v1[i] -= v2[i];
	}
	return v1;
}

vector<double>& operator-=(vector<double>& v1, double d)
{
	for (unsigned int i = 0; i < v1.size(); i++)
	{
		v1[i] -= d;
	}
	return v1;
}

vector<double>& operator*=(vector<double>& v, double x)
{
	for (unsigned int i = 0; i < v.size(); i++)
	{
		v[i] *= x;
	}
	return v;
}

vector<double>& operator*=(vector<double>& v, const vector<double>& v2)
{
	for (unsigned int i = 0; i < v.size(); i++)
	{
		v[i] *= v2[i];
	}
	return v;
}

vector<double>& operator/=(vector<double>& v, double x)
{
	for (unsigned int i = 0; i < v.size(); i++)
	{
		v[i] /= x;
	}
	return v;
}

vector<double> operator+(vector<double> v1, const vector<double>& v2)
{
	v1 += v2;
	return v1;
}

vector<double> operator+(vector<double> v1, double d)
{
	v1 += d;
	return v1;
}

vector<double> operator-(vector<double> v1, const vector<double>& v2)
{
	v1 -= v2;
	return v1;
}

vector<double> operator-(vector<double> v1, double d)
{
	v1 -= d;
	return v1;
}

vector<double> operator*(vector<double> v, double x)
{
	v *= x;
	return v;
}

vector<double> operator/(vector<double> v, double x)
{
	v /= x;
	return v;
}

vector<vector<double> >& operator+=(vector<vector<double> >& v1, const vector<vector<double> >& v2)
{
	for (unsigned int i = 0; i < v1.size(); i++)
	{
		v1[i] += v2[i];
	}
	return v1;
}

vector<vector<double> >& operator*=(vector<vector<double> >& v, double x)
{
	for (unsigned int i = 0; i < v.size(); i++)
	{
		v[i] *= x;
	}
	return v;
}

vector<vector<double> >& operator/=(vector<vector<double> >& v, double x)
{
	for (unsigned int i = 0; i < v.size(); i++)
	{
		v[i] /= x;
	}
	return v;
}

vector<vector<double> > operator+(vector<vector<double> > v1, const vector<vector<double> >& v2)
{
	v1 += v2;
	return v1;
}

vector<vector<double> > operator*(vector<vector<double> > v, double x)
{
	v *= x;
	return v;
}

vector<vector<double> > operator/(vector<vector<double> > v, double x)
{
	v /= x;
	return v;
}


Observables::ObservableCollection operator+(Observables::ObservableCollection lhs, const Observables::ObservableCollection& rhs)
{
	auto result = lhs.Clone();
	result += rhs;
	return result;
}

Observables::ObservableCollection operator-(Observables::ObservableCollection lhs, const Observables::ObservableCollection& rhs)
{
	auto result = lhs.Clone();
	result -= rhs;
	return result;
}

Observables::ObservableCollection operator*(Observables::ObservableCollection lhs, double rhs)
{
	auto result = lhs.Clone();
	result *= rhs;
	return result;
}

Observables::ObservableCollection operator*(Observables::ObservableCollection lhs, Observables::ObservableCollection& rhs)
{
	auto result = lhs.Clone();
	result *= rhs;
	return result;
}

Observables::ObservableCollection operator/(Observables::ObservableCollection lhs, double rhs)
{
	auto result = lhs.Clone();
	result /= rhs;
	return result;
}
