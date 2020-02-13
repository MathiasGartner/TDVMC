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

class BosonsDiscrete: public IPhysicalSystem
{
private:
	int numberOfSplines;
	vector<double> splineSums; //indices: k (bin); for <O_k>
	int orderFD;
	vector<double> factorsFD;
	vector<vector<vector<double> > > splineSumsFDF;
	vector<vector<vector<double> > > splineSumsFDB;

	double halfLength;
	double maxDistance;

	int numOfOtherLocalOperators;
	vector<double> otherLocalOperators;

	vector<vector<vector<double> > > kValues;
	vector<double> kNorms;
	int numOfkValues;

	vector<double> splineSumsNew; //indices: k (bin); for <O_k>
	vector<double> sumOldPerBin;
	vector<double> sumNewPerBin;
	int changedParticleIndex;

	vector<double> nodes;

	Observables::ObservableVsOnGridWithScaling pairDistribution;
	Observables::ObservableVsOnGrid structureFactor;

private:
	double GetExternalPotential(vector<double>& r);

	vector<double> CalculateWFChange(vector<vector<double> >& R, int changedParticleIndex, int coord, int steps);

public:
	void SetNodes(vector<double> n);

//Implementation of IPhysicalSystem
public:
	BosonsDiscrete(vector<double>& params, string configDirectory);

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
