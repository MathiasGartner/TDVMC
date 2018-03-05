#pragma once

#include "Constants.h"
#include <iostream>
#include "IPhysicalSystem.h"
#include "MathOperators.h"
#include "Utils.h"

using namespace std;

class HeBulk: public IPhysicalSystem
{
private:
	string configDirectory;

	int	numberOfSplines;
	vector<double> splineSums; //indices: k (bin); for <O_k>
	vector<vector<vector<double> >  > splineSumsD; //indices: k (bin), n (particle), a (coordinate)
	vector<vector<double> > splineSumsD2; //indices: k (bin), n (particle)
	double mcMillanSum;
	vector<vector<double> > mcMillanSumD;
	vector<double> mcMillanSumD2;
	double rijSplit;
	double halfLength;
	double maxDistance;
	double nodePointSpacing;
	double nodePointSpacing2;

	int grBinCount;
	int grBinStartIndex;
	vector<double> grBins;
	vector<double> grBinVolumes;
	double grMaxDistance;
	double grNodePointSpacing;
	vector<vector<vector<double> > > kValues;
	int numOfkValues;

	double exponentNew;
	double mcMillanSumNew;
	vector<double> splineSumsNew; //indices: k (bin); for <O_k>
	vector<double> sumOldPerBin;
	vector<double> sumNewPerBin;

	double factorFirstSpline1;
	double factorFirstSpline2;
	double factorSecondSpline1;
	double factorSecondSpline2;
	double factorLastSpline;
	double factorSecondLastSpline;
	double factorLastSplinePhi;
	double factorSecondLastSplinePhi;

//Implementation of IPhysicalSystem
public:
	HeBulk(string configDirectory);

	void InitSystem();

	vector<double> GetCenterOfMass(double** R);

	void CalculateExpectationValues(double** R, double* uR, double* uI, double phiR, double phiI);

	void CalculateAdditionalSystemProperties(double** R, double* uR, double* uI, double phiR, double phiI);

	void CalculateWavefunction(double** R, double* uR, double* uI, double phiR, double phiI);

	void CalculateWFChange(double** R, double* uR, double* uI, double phiR, double phiI, int changedParticleIndex, double* oldPosition);

	double GetWFQuotient(double** R, double* uR, double* uI, double phiR, double phiI, int changedParticleIndex, double* oldPosition);

	void AcceptMove();
};

