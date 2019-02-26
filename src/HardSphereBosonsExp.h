#pragma once

#include "CSDataBulkSplines.h"
#include "Constants.h"
#include "IPhysicalSystem.h"
#include "MathOperators.h"
#include "Utils.h"

#include <iostream>

using namespace std;

class HardSphereBosonsExp: public IPhysicalSystem
{
private:
	double halfLength;
	double maxDistance;

	int numOfOtherLocalOperators;
	vector<double> otherLocalOperators;
	int grBinCount;
	int grBinStartIndex;
	vector<double> grBins;
	vector<double> grBinVolumes;
	double grMaxDistance;
	double grNodePointSpacing;
	vector<vector<vector<double> > > kValues;
	int numOfkValues;

	vector<double> splineSumsNew; //indices: k (bin); for <O_k>
	vector<double> sumOldPerBin;
	vector<double> sumNewPerBin;
	int changedParticleIndex;

	vector<vector<double> > bcFactors; //factors according to the boundary conditions

protected:
	double exponentNew;
	vector<double> scalingFactors;
	int Mode;

//Implementation of IPhysicalSystem
public:
	HardSphereBosonsExp(vector<double>& params, string configDirectory);

	void InitSystem();

	vector<double> GetCenterOfMass(vector<vector<double> >& R);

	void CalculateOtherLocalOperators(vector<vector<double> >& R);

	void CalculateExpectationValues(vector<vector<double> >& R, vector<double>& uR, vector<double>& uI, double phiR, double phiI);

	void CalculateExpectationValues(ICorrelatedSamplingData* sample, vector<double>& uR, vector<double>& uI, double phiR, double phiI);

	void CalculateAdditionalSystemProperties(vector<vector<double> >& R, vector<double>& uR, vector<double>& uI, double phiR, double phiI);

	void CalculateWavefunction(vector<vector<double> >& R, vector<double>& uR, vector<double>& uI, double phiR, double phiI);

	void CalculateWavefunction(ICorrelatedSamplingData* sample, vector<double>& uR, vector<double>& uI, double phiR, double phiI);

	void CalculateWFChange(vector<vector<double> >& R, vector<double>& uR, vector<double>& uI, double phiR, double phiI, int changedParticleIndex, vector<double>& oldPosition);

	double CalculateWFQuotient(vector<vector<double> >& R, vector<double>& uR, vector<double>& uI, double phiR, double phiI, int changedParticleIndex, vector<double>& oldPosition);

	void AcceptMove();

	void FillCorrelatedSamplingData(ICorrelatedSamplingData* data);
};
