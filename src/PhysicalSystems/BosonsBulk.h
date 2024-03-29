#pragma once

#include "IPhysicalSystem.h"

#include "../Observables/ObservableVsOnGridWithScaling.h"

#include "../Constants.h"
#include "../CSDataBulkSplines.h"
#include "../MathOperators.h"
#include "../Utils.h"

#include <iostream>

using namespace std;

namespace PhysicalSystems
{

class BosonsBulk: public IPhysicalSystem
{
private:
	int numberOfSplines;
	double outerSum;
	vector<double> splineSums; //indices: k (bin); for <O_k>
	vector<vector<vector<double> > > splineSumsD; //indices: k (bin), n (particle), a (coordinate)
	vector<vector<double> > splineSumsD2; //indices: k (bin), n (particle)

	double halfLength;
	double maxDistance;
	int numberOfSpecialParametersStart;
	int numberOfSpecialParametersEnd;
	int numberOfStandardParameters;

	int numOfOtherLocalOperators;
	vector<double> otherLocalOperators;
	
	int numOfPairDistributionValues;
	vector<vector<vector<double> > > kValues;
	vector<double> kNorms;
	int numOfkValues;

	double outerSumNew;
	vector<double> splineSumsNew; //indices: k (bin); for <O_k>
	vector<double> sumOldPerBin;
	vector<double> sumNewPerBin;

	int np1;
	int np2;
	int np3;
	vector<vector<double> > bcFactorsStart; //factors according to the boundary conditions at origin
	vector<vector<double> > bcFactorsEnd; //factors according to the boundary conditions at L/2

	vector<double> nodes;
	vector<vector<vector<double> > > splineWeights;

	Observables::ObservableVsOnGridWithScaling pairDistribution;
	Observables::ObservableVsOnGrid structureFactor;
	Observables::ObservableVsOnGrid structureFactorCos;
	Observables::ObservableVsOnGrid structureFactorSin;

protected:
	vector<double> scalingFactors;

private:
	double GetExternalPotential(vector<double>& r);
	void RefreshLocalOperators();
	void CalculateLocalOperators(vector<vector<double> >& R);
	void CalculateExpectationValues(vector<double>& O, vector<vector<vector<double> > >& sD, vector<vector<double> >& sD2, vector<double>& otherO, vector<double>& gr, vector<double>& uR, vector<double>& uI, double phiR, double phiI);
	void CalculateWavefunction(vector<double>& O, vector<double>& uR, vector<double>& uI, double phiR, double phiI);

public:
	void SetNodes(vector<double> n);
	void SetPairDistributionBinCount(double n);

//Implementation of IPhysicalSystem
public:
	BosonsBulk(vector<double>& params, string configDirectory);

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
