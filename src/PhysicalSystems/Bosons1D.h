#pragma once

#include "IPhysicalSystem.h"

#include "../Observables/ObservableVsOnGridWithScaling.h"

#include "../Constants.h"
#include "../CSDataBulkSplines.h"
#include "../CSDataBosons1D.h"
#include "../MathOperators.h"
#include "../Utils.h"

#include <iostream>

using namespace std;

namespace PhysicalSystems
{

class Bosons1D: public IPhysicalSystem
{
private:
	int numberOfSplines;
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

	vector<double> splineSumsNew; //indices: k (bin); for <O_k>
	vector<double> sumOldPerBin;
	vector<double> sumNewPerBin;
	int changedParticleIndex; //TODO: is this needed?

	int np1;
	int np2;
	int np3;
	vector<vector<double> > bcFactorsStart; //factors according to the boundary conditions at origin
	vector<vector<double> > bcFactorsEnd; //factors according to the boundary conditions at L/2

	vector<double> nodes;
	vector<vector<vector<double> > > splineWeights;

	vector<double> nodeDiffs;
	vector<double> binSizePerNode;
	vector<double> binSizePerNode_R;
	vector<vector<vector<double> > > preCalcSplineValues;
	//double preCalcSplineValues[103][4][100000];
	bool usePreCalcSplineValues;

	Observables::Observable randNum;
	Observables::ObservableVsOnGridWithScaling pairDistribution;
	Observables::ObservableVsOnGrid structureFactor;
	Observables::ObservableVsOnGrid structureFactorCos;
	Observables::ObservableVsOnGrid structureFactorSin;

	Observables::ObservableV overlapToInitialState; //real part, imaginary part, absolute value

	vector<double> tmpLocalOperators;
	vector<double> tmpSplineSums;
	int particleStartIndexForReducedSampling;

	double exponentI;
	double exponentNewI;
	vector<double> initialParamsR;
	vector<double> initialParamsI;
	double initialParamPhiR;
	double initialParamPhiI;

	vector<CSDataBosons1D*> initialSamples;
	int sampleNo;
	bool calculateOverlap;

private:
	double GetExternalPotential(vector<double>& r);
	void RefreshLocalOperators();
	void RefreshLocalOperators(vector<double>& locO, vector<double>& sSums);
	void CalculateLocalOperators(vector<vector<double> >& R);
	void CalculateLocalOperators(vector<vector<double> >& R, vector<double>& locO, vector<double>& sSums);
	void CalculateExpectationValues(vector<double>& O, vector<vector<vector<double> > >& sD, vector<vector<double> >& sD2, vector<double>& otherO, vector<double>& gr, vector<double>& uR, vector<double>& uI, double phiR, double phiI);
	void CalculateWavefunction(vector<double>& O, vector<double>& uR, vector<double>& uI, double phiR, double phiI);
	void CalculateWavefunction(vector<vector<double> >& R, vector<double>& locO, vector<double>& sSums, vector<double>& uR, vector<double>& uI, double phiR, double phiI);

public:
	void SetNodes(vector<double> n);
	void SetPairDistributionBinCount(double n);
	void SetInitialParameters(vector<double>& uR, vector<double>& uI, double phiR, double phiI);
	void SetParticleStartIndexForReducedSampling(int n);

	double CalculateOBDMKernel(vector<double>& r, vector<vector<double> >& R, vector<double>& uR, vector<double>& uI, double phiR, double phiI);

	void AddToInitialSamples(vector<vector<double> >& R, Observables::ObservableCollection& data);
	void SetSampleNo(int sampleNo, bool calculateOverlap);

//Implementation of IPhysicalSystem
public:
	Bosons1D(vector<double>& params, string configDirectory);

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
