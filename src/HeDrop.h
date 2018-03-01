#pragma once

#include "Constants.h"
#include <iostream>
#include "IPhysicalSystem.h"
#include "MathOperators.h"
#include "Utils.h"

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

	int grBinCount = 200;
	int grBinStartIndex = 3;
	vector<double> grBins;
	vector<double> grBinVolumes;
	double grMaxDistance;
	double grNodePointSpacing;
	vector<vector<vector<double> > > kValues;
	int numOfkValues = 84;
	double densityProfileMaxDistance;
	int numOfDensityProfileValues = 200;

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
	bool USE_NORMALIZATION_AND_PHASE() = 0;

	void InitSystem() = 0;

	vector<double> GetCenterOfMass(double** R) = 0;

	void CalculateExpectationValues(double** R, double* uR, double* uI, double phiR, double phiI) = 0;

	void CalculateAdditionalSystemProperties(double** R, double* uR, double* uI, double phiR, double phiI) = 0;

	void CalculateWavefunction(double** R, double* uR, double* uI, double phiR, double phiI) = 0;

	void CalculateWFChange(double** R, double* uR, double* uI, double phiR, double phiI, int changedParticleIndex, double* oldPosition) = 0;

	double GetWFQuotient(double** R, double* uR, double* uI, double phiR, double phiI, int changedParticleIndex, double* oldPosition) = 0;

	void AcceptMove() = 0;
};

