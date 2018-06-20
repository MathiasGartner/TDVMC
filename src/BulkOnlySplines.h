#pragma once

#include "Constants.h"
#include "IPhysicalSystem.h"
#include "MathOperators.h"
#include "Utils.h"

#include <iostream>

using namespace std;

class BulkOnlySplines: public IPhysicalSystem
{
public:
	int numberOfSplines;
	vector<double> splineSums; //indices: k (bin); for <O_k>
	vector<vector<vector<double> > > splineSumsD; //indices: k (bin), n (particle), a (coordinate)
	vector<vector<double> > splineSumsD2; //indices: k (bin), n (particle)

	double halfLength;
	double maxDistance;
	double nodePointSpacing;
	double nodePointSpacing2;
	int numberOfSpecialParameters;
	int numberOfStandardParameters;

	int grBinCount;
	int grBinStartIndex;
	vector<double> grBins;
	vector<double> grBinVolumes;
	double grMaxDistance;
	double grNodePointSpacing;
	vector<vector<vector<double> > > kValues;
	int numOfkValues;

	double exponentNew;
	vector<double> splineSumsNew; //indices: k (bin); for <O_k>
	vector<double> sumOldPerBin;
	vector<double> sumNewPerBin;
	vector<vector<vector<double> > > sumPerBinPerParticle; //sumPerBinPerParticle[p1][p2][bin] -> p1...particle1, p2...particle2, bin...bin
	vector<vector<double> > sumNewPerBinForChangedParticle; //sumPerBinPerParticle[p1][p2][bin] -> p...particle, bin...bin
	int changedParticleIndex;

	vector<vector<double> > bcFactors; //factors according to the boundary conditions

public:
	vector<vector<vector<double> > > allLocalKineticEnergiesD1;
	vector<vector<double> > allLocalKineticEnergiesD2;
	vector<double> allER1;
	vector<double> allER2;
	vector<double> allER1new;
	vector<double> allER2new;

	vector<vector<double> > splineSumsOnlyD2;

//Implementation of IPhysicalSystem
public:
	BulkOnlySplines(vector<double>& params, string configDirectory);

	void InitSystem();

	vector<double> GetCenterOfMass(vector<vector<double> >& R);

	void CalculateExpectationValues(vector<vector<double> >& R, vector<double>& uR, vector<double>& uI, double phiR, double phiI);

	void CalculateAdditionalSystemProperties(vector<vector<double> >& R, vector<double>& uR, vector<double>& uI, double phiR, double phiI);

	void CalculateWavefunction(vector<vector<double> >& R, vector<double>& uR, vector<double>& uI, double phiR, double phiI);

	void CalculateWFChange(vector<vector<double> >& R, vector<double>& uR, vector<double>& uI, double phiR, double phiI, int changedParticleIndex, vector<double>& oldPosition);
	void CalculateWFChange2(vector<vector<double> >& R, vector<double>& uR, vector<double>& uI, double phiR, double phiI, int changedParticleIndex, vector<double>& oldPosition);

	double CalculateWFQuotient(vector<vector<double> >& R, vector<double>& uR, vector<double>& uI, double phiR, double phiI, int changedParticleIndex, vector<double>& oldPosition);

	void AcceptMove();
};

