#include "GaussianWavepacket.h"

using namespace std;

namespace PhysicalSystems
{

GaussianWavepacket::GaussianWavepacket(vector<double>& params, string configDirectory) :
		IPhysicalSystem(params, configDirectory)
{
	this->USE_NORMALIZATION_AND_PHASE = true;
	this->USE_NIC = false;
	this->USE_MOVE_COM_TO_ZERO = false;
}

void GaussianWavepacket::InitSystem()
{
	numOfOtherExpectationValues = 0;
	numOfAdditionalSystemProperties = 0;

	wf = 0.0;
	exponent = 0.0;

	localEnergyR = 0;
	localEnergyI = 0;
	localOperators.resize(N_PARAM);
	ClearVector(localOperators);
	localOperatorsMatrix.resize(N_PARAM);
	for (auto &row : localOperatorsMatrix)
	{
		row.resize(N_PARAM);
	}
	ClearVector(localOperatorsMatrix);
	localOperatorlocalEnergyR.resize(N_PARAM);
	ClearVector(localOperatorlocalEnergyR);
	localOperatorlocalEnergyI.resize(N_PARAM);
	ClearVector(localOperatorlocalEnergyI);
	otherExpectationValues.resize(numOfOtherExpectationValues);
	ClearVector(otherExpectationValues);
	additionalSystemProperties.resize(numOfAdditionalSystemProperties);
	ClearVector(additionalSystemProperties);

	wfNew = 0.0;
}

vector<double> GaussianWavepacket::GetCenterOfMass(vector<vector<double> >& R)
{
	vector<double> com = { 0.0, 0.0, 0.0 };
	return com;
}

void GaussianWavepacket::CalculateOtherLocalOperators(vector<vector<double> >& R)
{

}

void GaussianWavepacket::CalculateExpectationValues(vector<vector<double> >& R, vector<double>& uR, vector<double>& uI, double phiR, double phiI)
{
	double value = 0;
	double r = R[0][0];

	//E_R
	//kinetic
	value = pow(uR[0], 2) - pow(uI[0], 2) - 4 * r * uI[0] * uI[1] + 4 * r * uR[0] * uR[1] - 4 * pow(r, 2) * pow(uI[1], 2) + 4 * pow(r, 2) * pow(uR[1], 2) + 2 * uR[1];
	//for wavefunction with x and x^4
	//value = -pow(u[3],2) + pow(u[2],2) - 8*pow(r,3)*u[3]*u[5] + 8*pow(r,3)*u[2]*u[4] + 4*pow(r,2)*(-4*pow(r,4)*pow(u[5],2) + 3*u[4] + 4*pow(r,4)*pow(u[4],2));
	value *= -0.5; //-HBAR2_2M

	//potential
	double w = 1.0;
	value += (pow(w, 2) / 2.0) * pow(r, 2);
	//value += (pow(w, 2) / 2.0) * pow(r, 2) + (pow(w, 2) / 5.0) * pow(r, 4);

	localEnergyR = value;

	//E_I
	value = 2 * uR[0] * uI[0] + 4 * r * uR[0] * uI[1] + 4 * r * uI[0] * uR[1] + 8 * pow(r, 2) * uR[1] * uI[1] + 2 * uI[1];
	//for wavefunction with x and x^4
	//value = 2*(u[3]*(u[2] + 4*pow(r,3)*u[4]) + 2*pow(r,2)*u[5]*(3 + 2*r*u[2] + 8*pow(r,4)*u[4]));
	value *= -0.5; //-HBAR2_2M
	localEnergyI = value;

	localOperators[0] = r;
	localOperators[1] = r * r;

	for (int k = 0; k < N_PARAM; k++)
	{
		for (int j = 0; j < N_PARAM; j++)
		{
			localOperatorsMatrix[k][j] = localOperators[k] * localOperators[j];
		}
		localOperatorlocalEnergyR[k] = localOperators[k] * localEnergyR;
		localOperatorlocalEnergyI[k] = localOperators[k] * localEnergyI;
	}
}

void GaussianWavepacket::CalculateExpectationValues(ICorrelatedSamplingData* sample, vector<double>& uR, vector<double>& uI, double phiR, double phiI)
{

}

void GaussianWavepacket::CalculateAdditionalSystemProperties(vector<vector<double> >& R, vector<double>& uR, vector<double>& uI, double phiR, double phiI)
{

}

void GaussianWavepacket::CalculateWavefunction(vector<vector<double> >& R, vector<double>& uR, vector<double>& uI, double phiR, double phiI)
{
	//return psi^2
	double value = 1;
	value = exp(phiR) * exp(uR[0] * R[0][0]) * exp(uR[1] * pow(R[0][0], 2));
	//for wavefunction with x and x^4
	//value = exp(u[0]) * exp(u[2] * R[0][0]) * exp(u[4] * pow(R[0][0], 4));
	value = pow(value, 2);
	this->wf = value;
}

void GaussianWavepacket::CalculateWavefunction(ICorrelatedSamplingData* sample, vector<double>& uR, vector<double>& uI, double phiR, double phiI)
{

}

void GaussianWavepacket::CalculateWFChange(vector<vector<double> >& R, vector<double>& uR, vector<double>& uI, double phiR, double phiI, int changedParticleIndex, vector<double>& oldPosition)
{
	double tmp = this->wf;
	this->CalculateWavefunction(R, uR, uI, phiR, phiI);
	this->wfNew = this->wf;
	this->wf = tmp;
}

double GaussianWavepacket::CalculateWFQuotient(vector<vector<double> >& R, vector<double>& uR, vector<double>& uI, double phiR, double phiI, int changedParticleIndex, vector<double>& oldPosition)
{
	CalculateWFChange(R, uR, uI, phiR, phiI, changedParticleIndex, oldPosition);
	double wfQuotient = this->wfNew / this->wf;
	return wfQuotient;
}

void GaussianWavepacket::AcceptMove()
{
	this->wf = this->wfNew;
}

void GaussianWavepacket::FillCorrelatedSamplingData(ICorrelatedSamplingData* data)
{
}

}

