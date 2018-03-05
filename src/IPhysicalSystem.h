#pragma once

#include <vector>

using namespace std;

class IPhysicalSystem
{
protected:
	double time;
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
	bool USE_NORMALIZATION_AND_PHASE; //INFO: parameters phiR and phiI
	bool USE_NIC; //INFO: nearest image convention
	bool USE_MOVE_COM_TO_ZERO; //INFO: move center of mass to zero

	IPhysicalSystem()
	{
		this->time = 0;
	}

	virtual ~IPhysicalSystem()
	{
	}

	double GetTime() { return time; }
	void SetTime(double time) { this->time = time; }
	double GetWf() { return wf; }
	double GetWfNew() { return wfNew; }
	double GetExponent() { return exponent; }

	double GetLocalEnergyR() { return localEnergyR; }
	double GetLocalEnergyI() { return localEnergyI; }
	vector<double> GetLocalOperators() { return localOperators; }
	vector<vector<double> > GetLocalOperatorsMatrix() { return localOperatorsMatrix; }
	vector<double> GetLocalOperatorlocalEnergyR() { return localOperatorlocalEnergyR; }
	vector<double> GetLocalOperatorlocalEnergyI() { return localOperatorlocalEnergyI; }

	vector<double> GetOtherExpectationValues() { return otherExpectationValues; }
	vector<double> GetAdditionalSystemProperties() { return additionalSystemProperties; }
	int GetNumOfOtherExpectationValues() { return numOfOtherExpectationValues; }
	int GetNumOfAdditionalSystemProperties() { return numOfAdditionalSystemProperties; }

	virtual void InitSystem() = 0;

	virtual vector<double> GetCenterOfMass(double** R) = 0;

	virtual void CalculateExpectationValues(double** R, double* uR, double* uI, double phiR, double phiI) = 0;

	virtual void CalculateAdditionalSystemProperties(double** R, double* uR, double* uI, double phiR, double phiI) = 0;

	virtual void CalculateWavefunction(double** R, double* uR, double* uI, double phiR, double phiI) = 0;

	virtual void CalculateWFChange(double** R, double* uR, double* uI, double phiR, double phiI, int changedParticleIndex, double* oldPosition) = 0;

	virtual double GetWFQuotient(double** R, double* uR, double* uI, double phiR, double phiI, int changedParticleIndex, double* oldPosition) = 0;

	virtual void AcceptMove() = 0;
};
