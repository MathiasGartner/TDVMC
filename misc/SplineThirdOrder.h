#pragma once

#include "ISpline.h"

#include <vector>

using namespace std;

class SplineThirdOrder : public ISpline
{

//Implementation of ISpline
public:
	SplineThirdOrder() : ISpline(3)
{
}

	inline vector<double> value(vector<double> x);
	inline vector<double> d1(vector<double> x);
	inline vector<double> d2(vector<double> x);
};

vector<double> SplineThirdOrder::value(vector<double> x)
{
	return {
		1.0 / 6.0 * x[2],
		1.0 / 6.0 * (1.0 + 3.0 * x[0] + 3.0 * x[1] - 3.0 * x[2]),
		1.0 / 6.0 * (4.0 - 6.0 * x[1] + 3.0 * x[2]),
		-1.0 / 6.0 * (-1.0 + 3.0 * x[0] - 3.0 * x[1] + x[2])
	};
}

vector<double> SplineThirdOrder::d1(vector<double> x)
{
	return {

		1.0 / 2.0 * x[1],
		1.0 / 6.0 * (3.0 + 6.0 * x[0] - 9.0 * x[1]),
		1.0 / 6.0 * (-12.0 * x[0] + 9.0 * x[1]),
		-1.0 / 2.0 * (1.0 - 2.0 * x[0] + x[1])
	};
}

vector<double> SplineThirdOrder::d2(vector<double> x)
{
	return {
		x[0],
		1.0 / 6.0 * (6.0 - 18.0 * x[0]),
		1.0 / 6.0 * (-12.0 + 18.0 * x[0]),
		1.0 - x[0]
	};
}
