#pragma once

#include <vector>

using namespace std;

vector<double>& operator+=(vector<double>& v1, const vector<double>& v2);
vector<double>& operator-=(vector<double>& v1, const vector<double>& v2);
vector<double>& operator*=(vector<double>& v, double x);
vector<double>& operator/=(vector<double>& v, double x);
vector<double> operator+(vector<double> v1, const vector<double>& v2);
vector<double> operator-(vector<double> v1, const vector<double>& v2);
vector<double> operator*(vector<double> v, double x);
vector<double> operator/(vector<double> v, double x);

vector<vector<double> >& operator+=(vector<vector<double> >& v1, const vector<vector<double> >& v2);
vector<vector<double> >& operator*=(vector<vector<double> >& v, double x);
vector<vector<double> >& operator/=(vector<vector<double> >& v, double x);
vector<vector<double> > operator+(vector<vector<double> > v1, const vector<vector<double> >& v2);
vector<vector<double> > operator*(vector<vector<double> > v, double x);
vector<vector<double> > operator/(vector<vector<double> > v, double x);
