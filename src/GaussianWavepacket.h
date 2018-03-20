#pragma once

#include "Constants.h"
#include <iostream>
#include "IPhysicalSystem.h"
#include "MathOperators.h"
#include "Utils.h"

using namespace std;

class GaussianWavepacket: public IPhysicalSystem
{
private:
	string configDirectory;

//Implementation of IPhysicalSystem
public:
	GaussianWavepacket(string configDirectory);

	void InitSystem();

	vector<double> GetCenterOfMass(vector<vector<double> >& R);

	void CalculateExpectationValues(vector<vector<double> >& R, vector<double>& uR, vector<double>& uI, double phiR, double phiI);

	void CalculateAdditionalSystemProperties(vector<vector<double> >& R, vector<double>& uR, vector<double>& uI, double phiR, double phiI);

	void CalculateWavefunction(vector<vector<double> >& R, vector<double>& uR, vector<double>& uI, double phiR, double phiI);

	void CalculateWFChange(vector<vector<double> >& R, vector<double>& uR, vector<double>& uI, double phiR, double phiI, int changedParticleIndex, vector<double>& oldPosition);

	double CalculateWFQuotient(vector<vector<double> >& R, vector<double>& uR, vector<double>& uI, double phiR, double phiI, int changedParticleIndex, vector<double>& oldPosition);

	void AcceptMove();
};

