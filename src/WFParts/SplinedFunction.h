/*-----------------------------------------------------------------------------
 *
 * 		Name:			SplinedFunction.h
 * 		Author:			Mathias Gartner
 * 		Description:	Interface for a one-dimensional splined function
 * 						used for wavefunction modeling.
 * 						Holds the underlying grid, derivatives of the function,
 * 						as well as restictions due to boundary conditions on the
 * 						wavefunction.
 *
 * --------------------------------------------------------------------------*/

#pragma once

#include "../Utils.h"

using namespace std;

namespace WFParts
{

class SplinedFunction
{
public:
	int splineOrder;
	int numOfSplineParts;

	vector<double> nodes;
	double nodeSpacing;

	int numberOfSplines;
	vector<double> splineSums; //indices: k (bin); for <O_k>
	vector<vector<vector<double> > > splineSumsD; //indices: k (bin), n (particle), a (coordinate)
	vector<vector<double> > splineSumsD2; //indices: k (bin), n (particle)

	vector<double> splineSumsNew; //indices: k (bin); for <O_k>
	vector<double> sumOldPerBin;
	vector<double> sumNewPerBin;

	int numberOfSpecialParametersStart;
	int numberOfSpecialParametersEnd;
	int numberOfStandardParameters;
	int numberOfTotalParameters;
	int np1;
	int np2;
	int np3;
	vector<vector<double> > bcFactorsStart; //factors according to the boundary conditions at origin
	vector<vector<double> > bcFactorsEnd; //factors according to the boundary conditions at L/2

	vector<vector<vector<double> > > splineWeights; //indices: s (spline), p (part), c (coefficients for rij^n)

public:
	SplinedFunction(int splineOrder_);

	virtual void Init();
	int BinIndex(double r);
};

}
