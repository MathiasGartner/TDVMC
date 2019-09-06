#pragma once

#include "ICorrelatedSamplingData.h"

#include <vector>

using namespace std;

class CSDataBulkSplinesBR: public ICorrelatedSamplingData
{
public:
	vector<vector<vector<double> > > splineSumsD; //indices: k (bin), n (particle), a (coordinate)
	vector<vector<double> > splineSumsD2; //indices: k (bin), n (particle)
	vector<vector<vector<double> > > splineSumsDRad; //indices: k (bin), n (particle), a (coordinate)
	vector<vector<double> > splineSumsD2Rad; //indices: k (bin), n (particle)
	vector<double> otherLocalOperators; //INFO: potentialIntern, potentialInternComplex, potentialExtern, potentialExternComplex
	vector<double> grBins;
};
