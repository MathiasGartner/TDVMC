#pragma once

#include <vector>

#include "Utils.h"

using namespace std;

namespace SplineFactory
{

vector<vector<vector<double> > > GetWeights(vector<double> n);
void SetBoundaryConditions1_1(vector<double> n, vector<vector<double> >& bcs);
void SetBoundaryConditions2_2(vector<double> n, vector<vector<double> >& bcs);

//BC1.MM.1 Connection McMillan variable - Cubic NURBS
void SetBoundaryConditions1_MM_1(vector<double> n, vector<vector<double> >& bcs, double rc, double paramMcMillanExponent);
void SetBoundaryConditions1_EXP_1(vector<double> n, vector<vector<double> >& bcs, double rc);
void SetBoundaryConditions1_EXP_2(vector<double> n, vector<vector<double> >& bcs, double rc);
void SetBoundaryConditions1_EXP_3(vector<double> n, vector<vector<double> >& bcs, double rc);

}
