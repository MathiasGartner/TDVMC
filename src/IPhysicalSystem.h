#pragma once

#include <vector>

using namespace std;

class IPhysicalSystem
{
protected:

	vector<double> params;
	string configDirectory;

	double time;
	int step;
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

	IPhysicalSystem(vector<double>& params, string configDirectory)
	{
		this->params = params;
		this->configDirectory = configDirectory;

		this->time = 0;
		this->step = 0;
		this->wf = 0;
		this->wfNew = 0;
		this->exponent = 0;
		this->localEnergyR = 0;
		this->localEnergyI = 0;
		this->numOfOtherExpectationValues = 0;
		this->numOfAdditionalSystemProperties = 0;
		this->USE_NORMALIZATION_AND_PHASE = false;
		this->USE_NIC = false;
		this->USE_MOVE_COM_TO_ZERO = false;
	}

	virtual ~IPhysicalSystem()
	{
	}

	double GetTime() { return time; }
	void SetTime(double time) { this->time = time; }
	int GetStep() { return step; }
	void SetStep(double step) { this->step = step; }
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

	virtual vector<double> GetCenterOfMass(vector<vector<double> >& R) = 0;

	virtual void CalculateExpectationValues(vector<vector<double> >& R, vector<double>& uR, vector<double>& uI, double phiR, double phiI) = 0;

	virtual void CalculateAdditionalSystemProperties(vector<vector<double> >& R, vector<double>& uR, vector<double>& uI, double phiR, double phiI) = 0;

	virtual void CalculateWavefunction(vector<vector<double> >& R, vector<double>& uR, vector<double>& uI, double phiR, double phiI) = 0;

	virtual void CalculateWFChange(vector<vector<double> >& R, vector<double>& uR, vector<double>& uI, double phiR, double phiI, int changedParticleIndex, vector<double>& oldPosition) = 0;

	virtual double CalculateWFQuotient(vector<vector<double> >& R, vector<double>& uR, vector<double>& uI, double phiR, double phiI, int changedParticleIndex, vector<double>& oldPosition) = 0;

	virtual void AcceptMove() = 0;
};
