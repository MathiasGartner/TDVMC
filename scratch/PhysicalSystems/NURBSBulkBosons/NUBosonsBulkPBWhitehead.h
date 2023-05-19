#pragma once

#include "IPhysicalSystem.h"

#include "../Observables/ObservableVsOnGridWithScaling.h"
#include "../Observables/ObservableVsOnMultiGrid.h"

#include "../Constants.h"
#include "../CSDataBulkSplines.h"
#include "../MathOperators.h"
#include "../Utils.h"

#include <iostream>

using namespace std;

namespace PhysicalSystems
{

class NUBosonsBulkPBWhitehead: public IPhysicalSystem
{
private:
	int numberOfSplines;
	vector<double> splineSums; //indices: k (bin); for <O_k>
	vector<vector<vector<double> > > splineSumsD; //indices: k (bin), n (particle), a (coordinate)
	vector<vector<double> > splineSumsD2; //indices: k (bin), n (particle)

	double halfLength;
	double maxDistance;
	int numberOfSpecialParameters;
	int numberOfStandardParameters;

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

	vector<double> nodes;
	vector<vector<vector<double> > > splineWeights;

	Observables::ObservableVsOnGridWithScaling pairDistributionRad;
	Observables::ObservableVsOnMultiGrid pairDistribution;

private:
	double GetExternalPotential(vector<double>& r);
	void RefreshLocalOperators();
	void CalculateLocalOperators(vector<vector<double> >& R);
	void CalculateExpectationValues(vector<double>& O, vector<vector<vector<double> > >& sD, vector<vector<double> >& sD2, vector<double>& otherO, vector<double>& gr, vector<double>& uR, vector<double>& uI, double phiR, double phiI);
	void CalculateWavefunction(vector<double>& O, vector<double>& uR, vector<double>& uI, double phiR, double phiI);

	double TransformCoordinates(vector<double>& original, vector<double>& transformed);

public:
	void SetNodes(vector<double> n);
	void SetGrBinCount(double n);

	//vector<vector<double> > test;

//Implementation of IPhysicalSystem
public:
	NUBosonsBulkPBWhitehead(vector<double>& params, string configDirectory);

	void InitSystem() override;

	vector<double> GetCenterOfMass(vector<vector<double> >& R) override;

	void CalculateOtherLocalOperators(vector<vector<double> >& R) override;

	void CalculateExpectationValues(vector<vector<double> >& R, vector<double>& uR, vector<double>& uI, double phiR, double phiI) override;

	void CalculateExpectationValues(ICorrelatedSamplingData* sample, vector<double>& uR, vector<double>& uI, double phiR, double phiI) override;

	void CalculateAdditionalSystemProperties(vector<vector<double> >& R, vector<double>& uR, vector<double>& uI, double phiR, double phiI) override;

	void CalculateWavefunction(vector<vector<double> >& R, vector<double>& uR, vector<double>& uI, double phiR, double phiI) override;

	void CalculateWavefunction(ICorrelatedSamplingData* sample, vector<double>& uR, vector<double>& uI, double phiR, double phiI) override;

	void CalculateWFChange(vector<vector<double> >& R, vector<double>& uR, vector<double>& uI, double phiR, double phiI, int changedParticleIndex, vector<double>& oldPosition) override;

	double CalculateWFQuotient(vector<vector<double> >& R, vector<double>& uR, vector<double>& uI, double phiR, double phiI, int changedParticleIndex, vector<double>& oldPosition) override;

	void AcceptMove() override;

	void InitCorrelatedSamplingData(vector<ICorrelatedSamplingData*>& data) override;

	void FillCorrelatedSamplingData(ICorrelatedSamplingData* data) override;
};

}
