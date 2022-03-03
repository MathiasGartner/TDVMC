#pragma once

#include "IPhysicalSystem.h"

#include "../Observables/Observable.h"
#include "../Observables/ObservableVsOnGrid.h"
#include "../Observables/ObservableVsOnGridWithScaling.h"

#include "../Constants.h"
#include "../MathOperators.h"
#include "../OneParticleData.h"
#include "../ParticlePairProperties.h"
#include "../ParticleProperties.h"
#include "../Utils.h"
#include "../WFParts/PairCorrelation.h"

#include <iostream>

using namespace std;
using namespace WFParts;

namespace PhysicalSystems
{

class Bosons1DMixture: public IPhysicalSystem
{
protected:

	enum ParticleType
	{
		Boson = 0,
		Impurity,
		ParticleType_COUNT
	};

	vector<double> globalNodes;

	vector<int> particleTypeIndexMapping;
	vector<vector<int> > correlationIndexMapping;
	vector<ParticleType> originalParticleTypes;
	vector<int> particleTypes;
	vector<vector<int> > correlationTypes;

	vector<PairCorrelation> pcs;
	vector<ParticleProperties> particleProperties;
	vector<ParticlePairProperties> particlePairProperties;
	vector<OneParticleData> oneParticleData;

private:

	double halfLength;
	double maxDistance;

	int numOfOtherLocalOperators;
	vector<double> otherLocalOperators;

	int numOfPairDistributionValues;
	vector<vector<vector<double> > > kValues;
	vector<double> kNorms;
	int numOfkValues;

	Observables::ObservableVsOnGridWithScaling pairDistribution;
	Observables::ObservableVsOnGrid structureFactor;

private:
	double GetExternalPotential(vector<double>& r);
	void RefreshLocalOperators();
	void CalculateLocalOperators(vector<vector<double> >& R);
	void CalculateExpectationValues(vector<double>& O, vector<double>& otherO, vector<double>& gr, vector<double>& uR, vector<double>& uI, double phiR, double phiI);
	void CalculateWavefunction(vector<double>& O, vector<double>& uR, vector<double>& uI, double phiR, double phiI);

public:
	void ExtendNodes(vector<double>& n);
	void SetNodes(vector<double> n);
	void SetPairDistributionBinCount(double n);
	void SetParticleType(vector<int> p);
	void SetParticleType(vector<ParticleType> p);

	~Bosons1DMixture();

//Implementation of IPhysicalSystem
public:
	Bosons1DMixture(vector<double>& params, string configDirectory);

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
