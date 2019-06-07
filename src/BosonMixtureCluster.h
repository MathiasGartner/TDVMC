#pragma once

#include "Constants.h"
#include "CorrelationFunctionData.h"
#include "ICorrelatedSamplingData.h"
#include "IPhysicalSystem.h"
#include "MathOperators.h"
#include "OneParticleData.h"
#include "ParticlePairProperties.h"
#include "ParticleProperties.h"
#include "Utils.h"

#include <iostream>

using namespace std;

class BosonMixtureCluster: public IPhysicalSystem
{
protected:

	enum ParticleType
	{
		He3 = 0,
		He4,
		Na,
		Li,
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

protected:
	double exponentNew;

public:
	void SetNodes(vector<double> n); //TODO: do this for different CorrelationTypes
	void SetDensityProfileBinCount(double n); //TODO: do this for different CorrelationTypes
	void SetParticleType(vector<int> p);
	void SetParticleType(vector<ParticleType> p);

	~BosonMixtureCluster();

//Implementation of IPhysicalSystem
public:
	BosonMixtureCluster(vector<double>& params, string configDirectory);

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

