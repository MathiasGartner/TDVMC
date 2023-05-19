#include "HardSphereBosonsExp.h"

namespace PhysicalSystems
{

HardSphereBosonsExp::HardSphereBosonsExp(vector<double>& params, string configDirectory) :
		IPhysicalSystem(params, configDirectory)
{
	this->USE_NORMALIZATION_AND_PHASE = false;
	this->USE_NIC = true;
	this->USE_MOVE_COM_TO_ZERO = false;

	halfLength = 0;
	maxDistance = 0;

	numOfOtherLocalOperators = 0;
	grBinCount = 0;
	grBinStartIndex = 0;
	grMaxDistance = 0;
	grNodePointSpacing = 0;
	numOfkValues = 0;

	exponentNew = 0;
	changedParticleIndex = 0;

	Mode = 1;
}

void HardSphereBosonsExp::InitSystem()
{
	halfLength = LBOX / 2.0;
	maxDistance = halfLength;

	wf = 0.0;

	localEnergyR = 0;

	wfNew = 0.0;

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

	grBinCount = 100;
	grBinStartIndex = 3;
	grBins.resize(grBinCount);
	ClearVector(grBins);
	grMaxDistance = maxDistance;
	grNodePointSpacing = grMaxDistance / (double) grBinCount;

	grBinVolumes.resize(grBinCount);
	ClearVector(grBinVolumes);
	for (int i = 0; i < grBinCount; i++)
	{
		grBinVolumes[i] = 4.0 * M_PI * pow(grNodePointSpacing * (i + 1), 3.0) / 3.0;
	}
	for (int i = grBinCount - 1; i > 0; i--)
	{
		grBinVolumes[i] = grBinVolumes[i] - grBinVolumes[i - 1];
	}


	numOfOtherExpectationValues = 3 + grBinCount;
	otherExpectationValues.resize(numOfOtherExpectationValues);
	ClearVector(otherExpectationValues);
}

void HardSphereBosonsExp::CalculateOtherLocalOperators(vector<vector<double> >& R)
{
	//INFO: not needed
}

vector<double> HardSphereBosonsExp::GetCenterOfMass(vector<vector<double> >& R)
{
	vector<double> com = { 0.0, 0.0, 0.0 };
	return com;
}

void HardSphereBosonsExp::CalculateExpectationValues(vector<vector<double> >& R, vector<double>& uR, vector<double>& uI, double phiR, double phiI)
{
	double grBinInterval;
	int grBin;
	double kineticR = 0;
	double out = 0;
	vector<double> term1;
	term1.resize(3);
	ClearVector(term1);
	double term2 = 0;
	double term3 = 0;
	ClearVector(grBins);
	ClearVector(otherExpectationValues);
	for (int k = 0; k < N; k++)
	{
		ClearVector(term1);
		term2 = 0;
		term3 = 0;
		for (int j = 0; j < N; j++)
		{
			vector<double> vecr(DIM);
			double rkj = VectorDisplacementNIC(R[k], R[j], vecr);
			if (j != k)
			{
				if(rkj < maxDistance){
					double efct = exp(-(rkj-uR[0])/uR[1]);
					double efctone = 1 - efct;
					// term1 += (1/alpha)*(exp(-(boson.distanceNIC(k,j).norm()-a)/alpha)/(1-exp(-(boson.distanceNIC(k,j).norm()-a)/alpha)))*(boson.distanceNIC(k,j)/(boson.distanceNIC(k,j).norm()));
					// term2 += (1/(alpha*alpha))*(exp(-(boson.distanceNIC(k,j).norm()-a)/alpha)/(1-exp(-(boson.distanceNIC(k,j).norm()-a)/alpha))) + (exp(-2*(boson.distanceNIC(k,j).norm()-a)/alpha)/pow(1-exp(-(boson.distanceNIC(k,j).norm()-a)/alpha),2));
					// term3 += (1/alpha)*(2/boson.distanceNIC(k,j).norm())*(exp(-(boson.distanceNIC(k,j).norm()-a)/alpha)/(1-exp(-(boson.distanceNIC(k,j).norm()-a)/alpha)));
					term1 += vecr * (1/uR[1])*(1/efctone)*efct*(1/rkj);
					term2 += (1/(uR[1]*uR[1]))*((1/(efctone*efctone))*efct*efct + (efct/efctone));
					term3 += (2/(rkj*uR[1]*efctone))*efct;
				}
			}

			if (j < k && rkj < grMaxDistance)
			{
				grBinInterval = rkj / grNodePointSpacing;
				grBin = int(grBinInterval);
				grBins[grBin] += 1.0 / (double)pow(grBin + 1, 2);
				//grBins[grBin] += 1.0 / grBinVolumes[grBin];
			}
		}
		out += -VectorNorm2(term1) + term2 - term3;
	}
	kineticR = out;
	localEnergyR = kineticR;

	otherExpectationValues[0] = kineticR;
	otherExpectationValues[1] = 0;
	otherExpectationValues[2] = wf;
	for (int i = 0; i < grBinCount; i++)
	{
		otherExpectationValues[i + grBinStartIndex] = grBins[i];
	}
}

void HardSphereBosonsExp::CalculateExpectationValues(ICorrelatedSamplingData* sample, vector<double>& uR, vector<double>& uI, double phiR, double phiI)
{
	//INFO: not needed
}

void HardSphereBosonsExp::CalculateAdditionalSystemProperties(vector<vector<double> >& R, vector<double>& uR, vector<double>& uI, double phiR, double phiI)
{
	//INFO: not needed
}

void HardSphereBosonsExp::CalculateWavefunction(vector<vector<double> >& R, vector<double>& uR, vector<double>& uI, double phiR, double phiI)
{
	wf = 1;
}

void HardSphereBosonsExp::CalculateWavefunction(ICorrelatedSamplingData* sample, vector<double>& uR, vector<double>& uI, double phiR, double phiI)
{
	//INFO: not needed
}

void HardSphereBosonsExp::CalculateWFChange(vector<vector<double> >& R, vector<double>& uR, vector<double>& uI, double phiR, double phiI, int changedParticleIndex, vector<double>& oldPosition)
{
	double out = 1;
	double oldValues = 1;
	double jastrow = 0;
	for (int i = 0; i < N; i++)
	{
		if (i != changedParticleIndex)
		{
			if (Mode == 1)
			{
//				out*= 1+A;
				vector<double> vecrni(DIM);
				double r = VectorDisplacementNIC(R[i], R[changedParticleIndex], vecrni);
				if (r <= maxDistance)
				{
					if (r <= uR[0])
					{
						jastrow = 0;
					}
					else
					{
						jastrow = 1.0-exp(-(r - uR[0])/uR[1]);
					}
					out *= jastrow;
				}

				r = VectorDisplacementNIC(R[i], oldPosition, vecrni);
				if (r <= maxDistance)
				{
					if (r <= uR[0])
					{
						jastrow = 0;
					}
					else
					{
						jastrow = 1.0-exp(-(r - uR[0])/uR[1]);
					}
					oldValues *= jastrow;
				}
			}
		}
	}
	wfNew = out / oldValues;
}

double HardSphereBosonsExp::CalculateWFQuotient(vector<vector<double> >& R, vector<double>& uR, vector<double>& uI, double phiR, double phiI, int changedParticleIndex, vector<double>& oldPosition)
{
	CalculateWFChange(R, uR, uI, phiR, phiI, changedParticleIndex, oldPosition);
	//double wfQuotient = wfNew / wf;
	double wfQuotient = wfNew;

	return wfQuotient;
}

void HardSphereBosonsExp::AcceptMove()
{
	wf = wfNew;
}

void HardSphereBosonsExp::FillCorrelatedSamplingData(ICorrelatedSamplingData* data)
{
	//INFO: not needed
}

}
