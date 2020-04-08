#pragma once

#include "IPhysicalSystem.h"

#include "../Observables/Observable.h"
#include "../Observables/ObservableVsOnGrid.h"
#include "../Observables/ObservableVsOnGridWithScaling.h"

#include "../Constants.h"
#include "../CorrelationFunctionData.h"
#include "../MathOperators.h"
#include "../OneParticleData.h"
#include "../ParticlePairProperties.h"
#include "../ParticleProperties.h"
#include "../Utils.h"

#include <iostream>

using namespace std;

namespace PhysicalSystems
{

class BosonMixtureCluster: public IPhysicalSystem
{
protected:

	enum ParticleType
	{
		He3 = 0,
		He4,
		Na,
		Li,
		Cs,
		ParticleType_COUNT
	};

	double globalDensityProfileMaxDistance;
	int globalNumOfDensityProfileValues;
	vector<double> globalDensityProfileBins;
	double globalDensityProfileBinInterval;
	int globalDensityProfileBin;
	double globalDensityProfileNodePointSpacing;
	vector<double> globalDensityProfileBinVolumes;

	vector<double> globalNodes;

	vector<int> particleTypeIndexMapping;
	vector<vector<int> > correlationIndexMapping;
	vector<ParticleType> originalParticleTypes;
	vector<int> particleTypes;
	vector<vector<int> > correlationTypes;

	vector<CorrelationFunctionData> corrFuncData;
	vector<ParticleProperties> particleProperties;
	vector<ParticlePairProperties> particlePairProperties;
	vector<OneParticleData> oneParticleData;

	Observables::Observable r2;
	Observables::ObservableVsOnGrid angularDistribution;
	Observables::ObservableVsOnGridWithScaling densityFromCOM;
	Observables::ObservableVsOnGrid particleDistances;

public:
	void SetNodes(vector<double> n); //TODO: do this for different CorrelationTypes
	void SetDensityProfileBinCount(double n); //TODO: do this for different CorrelationTypes
	void SetParticleType(vector<int> p);
	void SetParticleType(vector<ParticleType> p);

	~BosonMixtureCluster();

//Implementation of IPhysicalSystem
public:
	BosonMixtureCluster(vector<double>& params, string configDirectory);

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
