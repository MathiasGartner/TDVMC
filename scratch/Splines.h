#pragma once

#include <vector>

using namespace std;

namespace Spliness
{

vector<double> value(vector<double> x)
{
	return {
		1.0 / 6.0 * x[2],
		1.0 / 6.0 * (1.0 + 3.0 * x[0] + 3.0 * x[1] - 3.0 * x[2]),
		1.0 / 6.0 * (4.0 - 6.0 * x[1] + 3.0 * x[2]),
		-1.0 / 6.0 * (-1.0 + 3.0 * x[0] - 3.0 * x[1] + x[2])
	};
}

vector<double> d1(vector<double> x)
{
	return {

		1.0 / 2.0 * x[1],
		1.0 / 6.0 * (3.0 + 6.0 * x[0] - 9.0 * x[1]),
		1.0 / 6.0 * (-12.0 * x[0] + 9.0 * x[1]),
		-1.0 / 2.0 * (1.0 - 2.0 * x[0] + x[1])
	};
}

vector<double> d2(vector<double> x)
{
	return {
		x[0],
		1.0 / 6.0 * (6.0 - 18.0 * x[0]),
		1.0 / 6.0 * (-12.0 + 18.0 * x[0]),
		1.0 - x[0]
	};
}

}
