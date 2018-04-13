#pragma once

#include "Constants.h"
#include "IPhysicalSystem.h"
#include "MathOperators.h"
#include "Utils.h"

#include <iostream>

using namespace std;

class BulkOnlySplinesQuadraticTail: public IPhysicalSystem
{
private:
	string configDirectory;

	int	numberOfSplines;
	vector<double> splineSums; //indices: k (bin); for <O_k>
	vector<vector<vector<double> >  > splineSumsD; //indices: k (bin), n (particle), a (coordinate)
	vector<vector<double> > splineSumsD2; //indices: k (bin), n (particle)
	double quadraticSum;
	vector<vector<double> > quadraticSumD;
	vector<double> quadraticSumD2;

	double halfLength;
	double maxDistance;
	double rijSplit;
	double nodePointSpacing;
	double nodePointSpacing2;
	int numberOfSpecialParameters;
	int numberOfStandardParameters;

	int grBinCount;
	int grBinStartIndex;
	vector<double> grBins;
	vector<double> grBinVolumes;
	double grMaxDistance;
	double grNodePointSpacing;
	vector<vector<vector<double> > > kValues;
	int numOfkValues;

	double exponentNew;
	double quadraticSumNew;
	vector<double> splineSumsNew; //indices: k (bin); for <O_k>
	vector<double> sumOldPerBin;
	vector<double> sumNewPerBin;

	double a13;
	double a14;
	double a15;
	double al;

//Implementation of IPhysicalSystem
public:
	BulkOnlySplinesQuadraticTail(string configDirectory);

	void InitSystem();

	vector<double> GetCenterOfMass(vector<vector<double> >& R);

	void CalculateExpectationValues(vector<vector<double> >& R, vector<double>& uR, vector<double>& uI, double phiR, double phiI);

	void CalculateAdditionalSystemProperties(vector<vector<double> >& R, vector<double>& uR, vector<double>& uI, double phiR, double phiI);

	void CalculateWavefunction(vector<vector<double> >& R, vector<double>& uR, vector<double>& uI, double phiR, double phiI);

	void CalculateWFChange(vector<vector<double> >& R, vector<double>& uR, vector<double>& uI, double phiR, double phiI, int changedParticleIndex, vector<double>& oldPosition);

	double CalculateWFQuotient(vector<vector<double> >& R, vector<double>& uR, vector<double>& uI, double phiR, double phiI, int changedParticleIndex, vector<double>& oldPosition);

	void AcceptMove();
};

