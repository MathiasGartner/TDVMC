#pragma once

#include "IPhysicalSystem.h"

#include "../Observables/ObservableVsOnGridWithScaling.h"

#include "../WFParts/PairCorrelation.h"
#include "../WFParts/SingleParticleFunction.h"

#include "../Constants.h"
#include "../CSDataBulkSplines.h"
#include "../MathOperators.h"
#include "../Utils.h"

#include <iostream>

using namespace std;
using namespace WFParts;

namespace PhysicalSystems
{

class InhBosons1D: public IPhysicalSystem
{
private:
	SingleParticleFunction spf;
	PairCorrelation pc;

	double halfLength;
	double maxDistance;

	int numOfOtherLocalOperators;
	vector<double> otherLocalOperators;

	vector<vector<vector<double> > > kValues;
	vector<double> kNorms;
	int numOfkValues;

	Observables::ObservableVsOnGridWithScaling density;
	Observables::ObservableVsOnGridWithScaling pairDistribution;
	Observables::ObservableVsOnGrid structureFactor;

private:
	double GetExternalPotential(vector<double>& r);
	void RefreshLocalOperators();
	void CalculateLocalOperators(vector<vector<double> >& R);
	void CalculateExpectationValues(vector<double>& O, vector<vector<vector<double> > >& spf_sD, vector<vector<double> >& spf_sD2, vector<vector<vector<double> > >& pc_sD, vector<vector<double> >& pc_sD2, vector<double>& otherO, vector<double>& uR, vector<double>& uI, double phiR, double phiI);
	void CalculateWavefunction(vector<double>& O, vector<double>& uR, vector<double>& uI, double phiR, double phiI);
	void ExtendNodes(vector<double>& n);

public:
	void SetNodesSPF(vector<double> n);
	void SetNodesPC(vector<double> n);

//Implementation of IPhysicalSystem
public:
	InhBosons1D(vector<double>& params, string configDirectory);

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
