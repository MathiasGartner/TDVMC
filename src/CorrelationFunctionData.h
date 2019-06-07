#pragma once

#include "Constants.h"
#include "Utils.h"

#include <vector>

using namespace std;

class CorrelationFunctionData
{
public:
	bool isUsed;

	int numberOfSplines;
	vector<double> splineSums; //indices: k (bin); for <O_k>
	vector<vector<vector<double> > > splineSumsD; //indices: k (bin), n (particle), a (coordinate)
	vector<vector<double> > splineSumsD2; //indices: k (bin), n (particle)
	double mcMillanSum;
	vector<vector<double> > mcMillanSumD;
	vector<double> mcMillanSumD2;
	double constSum;
	vector<vector<double> > constSumD;
	vector<double> constSumD2;
	double logSum;
	vector<vector<double> > logSumD;
	vector<double> logSumD2;
	double linearSum;
	vector<vector<double> > linearSumD;
	vector<double> linearSumD2;

	double rijSplit;
	double rijTail;
	double mcMillanFactor;

	double mcMillanSumNew;
	double constSumNew;
	double logSumNew;
	double linearSumNew;
	vector<double> splineSumsNew; //indices: k (bin); for <O_k>
	vector<double> sumOldPerBin;
	vector<double> sumNewPerBin;

	vector<vector<double> > bcFactors; //factors according to the boundary conditions

	vector<double> nodes;
	vector<vector<vector<double> > > splineWeights;

	//INFO: for WFChange calculation
	double mcMillanOld;
	double mcMillanNew;
	double constOld;
	double constNew;
	double logOld;
	double logNew;
	double linearOld;
	double linearNew;

public:
	CorrelationFunctionData();

	void Init();
};
