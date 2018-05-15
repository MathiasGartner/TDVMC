#pragma once

#include <vector>

using namespace std;

class ISpline
{
public:
	int order;

	ISpline(int order)
	{
		this->order = order;
	}

	virtual ~ISpline()
	{
	}

	virtual vector<double> value(vector<double> x) = 0;
	virtual vector<double> d1(vector<double> x) = 0;
	virtual vector<double> d2(vector<double> x) = 0;
};
