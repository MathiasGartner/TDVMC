#pragma once

#include <vector>

using namespace std;

class ICorrelatedSamplingData
{
public:
	vector<vector<double> > R;
	double wf;
	double wf2;
	double exponent;
	double exponent2;

	vector<double> localOperators; //indices: k (bin); for <O_k>

public:
	ICorrelatedSamplingData()
	{
		wf = 0;
		wf2 = 0;
		exponent = 0;
		exponent2 = 0;
	}

	virtual ~ICorrelatedSamplingData()
	{
	}
};
