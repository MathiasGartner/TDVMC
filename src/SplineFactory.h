#pragma once

#include <vector>

#include "Utils.h"

using namespace std;

namespace SplineFactory
{

vector<vector<vector<double> > > GetWeights(vector<double> n);
vector<vector<vector<double> > > GetWeights(vector<double> n, int splineOrder);
vector<vector<vector<double> > > GetWeights0(vector<double> n);
vector<vector<vector<double> > > GetWeights3(vector<double> n);
vector<vector<vector<double> > > GetWeights4(vector<double> n);

//for 1D systems with NURBS, 0th order
void SetBoundaryConditions0_1D_OR_1(vector<double> n, vector<vector<double> >& bcs);
void SetBoundaryConditions0_1D_CO_1(vector<double> n, vector<vector<double> >& bcs, double rc);

void SetBoundaryConditions1_1(vector<double> n, vector<vector<double> >& bcs);
void SetBoundaryConditions1_2(vector<double> n, vector<vector<double> >& bcs);
void SetBoundaryConditions2_1(vector<double> n, vector<vector<double> >& bcs);
void SetBoundaryConditions2_2(vector<double> n, vector<vector<double> >& bcs);

//BC1.MM.1 Connection McMillan variable - Cubic NURBS
void SetBoundaryConditions1_MM_1(vector<double> n, vector<vector<double> >& bcs, double rc, double paramMcMillanExponent);
void SetBoundaryConditions1_EXP_1(vector<double> n, vector<vector<double> >& bcs, double rc);
void SetBoundaryConditions1_EXP_2(vector<double> n, vector<vector<double> >& bcs, double rc);
void SetBoundaryConditions1_EXP_3(vector<double> n, vector<vector<double> >& bcs, double rc);

void SetBoundaryConditions1_MM_1_4thorder(vector<double> n, vector<vector<double> >& bcs, double rc, double paramMcMillanExponent);
void SetBoundaryConditions1_EXP_2_4thorder(vector<double> n, vector<vector<double> >& bcs, double rc);

//for 1D systems with NURBS, 3rd order
void SetBoundaryConditions3_1D_OR_1(vector<double> n, vector<vector<double> >& bcs, bool uniform = false);
void SetBoundaryConditions3_1D_OR_2(vector<double> n, vector<vector<double> >& bcs, bool uniform = false);
void SetBoundaryConditions3_1D_OR_3(vector<double> n, vector<vector<double> >& bcs, bool uniform = false);
void SetBoundaryConditions3_1D_OR_4(vector<double> n, vector<vector<double> >& bcs, bool uniform = false);
void SetBoundaryConditions3_1D_CO_1(vector<double> n, vector<vector<double> >& bcs, double rc, bool uniform = false);
void SetBoundaryConditions3_1D_CO_2(vector<double> n, vector<vector<double> >& bcs, double rc, bool uniform = false);
void SetBoundaryConditions3_1D_CO_3(vector<double> n, vector<vector<double> >& bcs, double rc, bool uniform = false);
void SetBoundaryConditions3_1D_CO_4(vector<double> n, vector<vector<double> >& bcs, double rc, bool uniform = false);

//for 1D systems with NURBS, 4th order
void SetBoundaryConditions4_1D_OR_1(vector<double> n, vector<vector<double> >& bcs);
void SetBoundaryConditions4_1D_CO_2(vector<double> n, vector<vector<double> >& bcs, double rc);
}
