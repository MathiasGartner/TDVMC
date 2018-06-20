#pragma once

#include "Constants.h"
#include "IPhysicalSystem.h"
#include "MathOperators.h"
#include "Utils.h"

#include <iostream>

using namespace std;

class HeDrop: public IPhysicalSystem
{
private:
	int numberOfSplines;
	int numberOfShortSplines;
	int numberOfLargeSplines;
	vector<double> splineSums; //indices: k (bin); for <O_k>
	vector<vector<vector<double> > > splineSumsD; //indices: k (bin), n (particle), a (coordinate)
	vector<vector<double> > splineSumsD2; //indices: k (bin), n (particle)
	double mcMillanSum;
	vector<vector<double> > mcMillanSumD;
	vector<double> mcMillanSumD2;
	double constSum;
	vector<vector<double> > constSumD;
	vector<double> constSumD2;
	//double logSum;
	//vector<vector<double> > logSumD;
	//vector<double> logSumD2;
	double linearSum;
	vector<vector<double> > linearSumD;
	vector<double> linearSumD2;
	double rijSplit;
	double rijSplineSplit;
	double nodePointSpacingShort;
	double nodePointSpacingShort2;
	double nodePointSpacingLarge;
	double nodePointSpacingLarge2;
	double rijTail;

	int grBinCount;
	int grBinStartIndex;
	vector<double> grBins;
	vector<double> grBinVolumes;
	double grMaxDistance;
	double grNodePointSpacing;
	vector<vector<vector<double> > > kValues;
	int numOfkValues;
	double densityProfileMaxDistance;
	int numOfDensityProfileValues;
	vector<double> densityProfileBins;
	double densityProfileBinInterval;
	int densityProfileBin;
	double densityProfileNodePointSpacing;

	double exponentNew;
	double mcMillanSumNew;
	double constSumNew;
	//double logSumNew;
	double linearSumNew;
	vector<double> splineSumsNew; //indices: k (bin); for <O_k>
	vector<double> sumOldPerBin;
	vector<double> sumNewPerBin;

	double factorFirstSpline1;
	double factorFirstSpline2;
	double factorSecondSpline1;
	double factorSecondSpline2;
	double factorSecondLastSpline;
	double factorSecondLastSplineConst;
	//double factorSecondLastSplineLog;
	double factorSecondLastSplineLinear;
	double factorLastSpline;
	double factorLastSplineConst;
	//double factorLastSplineLog;
	double factorLastSplineLinear;
	double factorSecondLastShort;
	double factorSecondLastLarge;
	double factorLastShort;
	double factorLastLarge;
	double factorFirstShort;
	double factorFirstLarge;
	double factorSecondShort;
	double factorSecondLarge;

//Implementation of IPhysicalSystem
public:
	HeDrop(vector<double>& params, string configDirectory);

	void InitSystem();

	vector<double> GetCenterOfMass(vector<vector<double> >& R);

	void CalculateExpectationValues(vector<vector<double> >& R, vector<double>& uR, vector<double>& uI, double phiR, double phiI);

	void CalculateAdditionalSystemProperties(vector<vector<double> >& R, vector<double>& uR, vector<double>& uI, double phiR, double phiI);

	void CalculateWavefunction(vector<vector<double> >& R, vector<double>& uR, vector<double>& uI, double phiR, double phiI);

	void CalculateWFChange(vector<vector<double> >& R, vector<double>& uR, vector<double>& uI, double phiR, double phiI, int changedParticleIndex, vector<double>& oldPosition);

	double CalculateWFQuotient(vector<vector<double> >& R, vector<double>& uR, vector<double>& uI, double phiR, double phiI, int changedParticleIndex, vector<double>& oldPosition);

	void AcceptMove();
};

