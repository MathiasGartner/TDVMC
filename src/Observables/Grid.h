#pragma once

#include <string>
#include <vector>

using namespace std;

namespace Observables
{

class Grid
{
public:
	string name;
	double min;
	double max;
	int count;
	double spacing;
	vector<double> grid;

public:
	Grid();

	void Init(double min, double max, double spacing);
	void Init(vector<double> g);

	int GetIndex(double value);
};

}
