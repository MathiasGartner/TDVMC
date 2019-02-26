#pragma once

#include <vector>

#include "Utils.h"

using namespace std;

namespace SplineFactory
{

vector<vector<vector<double> > > GetWeights(vector<double> n);
vector<double> GetBoundaryConditions1_1(vector<double> n);

}
