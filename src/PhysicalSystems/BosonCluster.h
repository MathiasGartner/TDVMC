#pragma once

#include "IPhysicalSystem.h"

#include "../Constants.h"
#include "../MathOperators.h"
#include "../Utils.h"

#include <iostream>

using namespace std;

namespace PhysicalSystems
{

class BosonCluster: public IPhysicalSystem
{
protected:
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
	//double logSum;
	//vector<vector<double> > logSumD;
	//vector<double> logSumD2;
	double linearSum;
	vector<vector<double> > linearSumD;
	vector<double> linearSumD2;

	double rijSplit;
	double rijTail;
	double mcMillanFactor;

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

	double mcMillanSumNew;
	double constSumNew;
	//double logSumNew;
	double linearSumNew;
	vector<double> splineSumsNew; //indices: k (bin); for <O_k>
	vector<double> sumOldPerBin;
	vector<double> sumNewPerBin;

	vector<vector<double> > bcFactors; //factors according to the boundary conditions

	vector<double> nodes;
	vector<vector<vector<double> > > splineWeights;

protected:
	double exponentNew;
	vector<double> scalingFactors;

public:
	void SetNodes(vector<double> n);
	void SetDensityProfileBinCount(double n);

//Implementation of IPhysicalSystem
public:
	BosonCluster(vector<double>& params, string configDirectory);

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

}
