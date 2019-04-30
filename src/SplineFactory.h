#pragma once

#include <vector>

#include "Utils.h"

using namespace std;

namespace SplineFactory
{

vector<vector<vector<double> > > GetWeights(vector<double> n);
void SetBoundaryConditions1_1(vector<double> n, vector<vector<double> >& bcs);
void SetBoundaryConditions2_2(vector<double> n, vector<vector<double> >& bcs);

}
