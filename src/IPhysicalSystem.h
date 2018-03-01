#pragma once

#include <vector>

using namespace std;

class IPhysicalSystem
{
protected:
	double wf;
	double wfNew;
	double exponent;

	double localEnergyR; // for <E^R>
	double localEnergyI; // for <E^I>
	vector<double> localOperators; //indices: k (bin); for <O_k>
	vector<vector<double> > localOperatorsMatrix; // for <O_k O_j>
	vector<double> localOperatorlocalEnergyR; // for <O_k E^R>
	vector<double> localOperatorlocalEnergyI; // for <O_k E^I>

	vector<double> otherExpectationValues; // eg. for potential and kinetic energy
	int numOfOtherExpectationValues;
	vector<double> additionalSystemProperties; // for properties at the end of the simulation
	int numOfAdditionalSystemProperties;

public:
	IPhysicalSystem()
	{
	}

	virtual ~IPhysicalSystem()
	{
	}

	virtual double GetWf() { return wf; }
	virtual double GetWfNew() { return wfNew; }
	virtual double GetExponent() { return exponent; }

	virtual double GetLocalEnergyR() { return localEnergyR; }
	virtual double GetLocalEnergyI() { return localEnergyI; }
	virtual vector<double> GetLocalOperators() { return localOperators; }
	virtual vector<vector<double> > GetLocalOperatorsMatrix() { return localOperatorsMatrix; }
	virtual vector<double> GetLocalOperatorlocalEnergyR() { return localOperatorlocalEnergyR; }
	virtual vector<double> GetLocalOperatorlocalEnergyI() { return localOperatorlocalEnergyI; }

	virtual vector<double> GetOtherExpectationValues() { return otherExpectationValues; }
	virtual vector<double> GetAdditionalSystemProperties() { return additionalSystemProperties; }
	virtual int GetNumOfOtherExpectationValues() { return numOfOtherExpectationValues; }
	virtual int GetNumOfAdditionalSystemProperties() { return numOfAdditionalSystemProperties; }

	virtual bool USE_NORMALIZATION_AND_PHASE() = 0;

	virtual void InitSystem() = 0;

	virtual vector<double> GetCenterOfMass(double** R) = 0;

	virtual void CalculateExpectationValues(double** R, double* uR, double* uI, double phiR, double phiI) = 0;

	virtual void CalculateAdditionalSystemProperties(double** R, double* uR, double* uI, double phiR, double phiI) = 0;

	virtual void CalculateWavefunction(double** R, double* uR, double* uI, double phiR, double phiI) = 0;

	virtual void CalculateWFChange(double** R, double* uR, double* uI, double phiR, double phiI, int changedParticleIndex, double* oldPosition) = 0;

	virtual double GetWFQuotient(double** R, double* uR, double* uI, double phiR, double phiI, int changedParticleIndex, double* oldPosition) = 0;

	virtual void AcceptMove() = 0;
};
