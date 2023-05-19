#include "BosonsDiscrete.h"

#include "../SplineFactory.h"

namespace PhysicalSystems
{

BosonsDiscrete::BosonsDiscrete(vector<double>& params, string configDirectory) :
		IPhysicalSystem(params, configDirectory)
{
	this->USE_NORMALIZATION_AND_PHASE = true;
	this->USE_NIC = true;
	this->USE_MOVE_COM_TO_ZERO = false;

	numberOfSplines = 0;

	halfLength = 0;
	maxDistance = 0;

	orderFD = 0;

	numOfOtherLocalOperators = 0;
	numOfkValues = 0;

	exponentNew = 0;
	changedParticleIndex = 0;
}

void BosonsDiscrete::SetNodes(vector<double> n)
{
	this->nodes = n;
	int periodicNodeCount = 0;
	for (int i = 0; i < periodicNodeCount; i++)
	{
		this->nodes.insert(this->nodes.begin(), -n[i + 1]);
		this->nodes.push_back(2.0 * n[n.size() - 1] - n[n.size() - 2 - i]);
	}
}

void BosonsDiscrete::InitSystem()
{
	numOfOtherLocalOperators = 4; //INFO: potentialIntern, potentialInternComplex, potentialExtern, potentialExternComplex
	otherLocalOperators.resize(numOfOtherLocalOperators);
	numOfkValues = 50;
	numOfOtherExpectationValues = 9;
	numOfAdditionalSystemProperties = numOfOtherExpectationValues;

	//INFO: BC
	//numberOfSplines = N_PARAM + 1; //BC3.1D.CO.1
	//numberOfSplines = N_PARAM + 2; //BC3.1D.CO.2
	//numberOfSplines = N_PARAM + 3; //BC3.1D.CO.2 + BC3.1D.OR.1
	//numberOfSplines = nodes.size() - 3 - 1; //3rd order splines -(d+1)
	numberOfSplines = nodes.size() - 1;
	halfLength = LBOX / 2.0;
	maxDistance = nodes[nodes.size() - 1];

	if (nodes.empty())
	{
		for (int i = -3; i < N_PARAM + 3; i++)
		{
			nodes.push_back((i * halfLength) / ((double) N_PARAM - 1));
		}
	}

	if (N_PARAM != numberOfSplines - 0)
	{
		cout << "!!! WRONG NUMBER OF PARAMETERS !!!" << endl;
	}

	wf = 0.0;
	exponent = 0.0;
	InitVector(splineSums, numberOfSplines, 0.0);
	//InitVector(splineSumsD, numberOfSplines, N, DIM, 0.0);
	//InitVector(splineSumsD2, numberOfSplines, N, 0.0);
	//orderFD = 5;
	//factorsFD = { -1.0/12.0, 4.0/3.0, -5.0/2.0, 4.0/3.0, -1.0/12.0 };
	orderFD = 1;
	factorsFD = { 1.0, -2.0, 1.0 };
	InitVector(splineSumsFDF, N, orderFD, numberOfSplines, 0.0);
	InitVector(splineSumsFDB, N, orderFD, numberOfSplines, 0.0);

	localEnergyR = 0;
	localEnergyI = 0;
	InitVector(localOperators, N_PARAM, 0.0);
	InitVector(localOperatorsMatrix, N_PARAM, N_PARAM, 0.0);
	InitVector(localOperatorlocalEnergyR, N_PARAM, 0.0);
	InitVector(localOperatorlocalEnergyI, N_PARAM, 0.0);
	InitVector(otherExpectationValues, numOfOtherExpectationValues, 0.0);
	InitVector(additionalSystemProperties, numOfAdditionalSystemProperties, 0.0);

	wfNew = 0.0;
	exponentNew = 0.0;
	InitVector(splineSumsNew, numberOfSplines, 0.0);
	InitVector(sumOldPerBin, numberOfSplines, 0.0);
	InitVector(sumNewPerBin, numberOfSplines, 0.0);

	kValues = ReadKValuesFromJsonFile(this->configDirectory + "kVectors" + to_string(DIM) + "D.json");
	ReadDataFromFile(kNorms, this->configDirectory + "kNorm" + to_string(DIM) + "D.csv");
	kNorms.resize(numOfkValues);
	for (int k = 0; k < numOfkValues; k++)
	{
		for (unsigned int kn = 0; kn < kValues[k].size(); kn++)
		{
			for (int a = 0; a < DIM; a++)
			{
				kValues[k][kn][a] *= 2 * M_PI / LBOX;
			}
		}
		kNorms[k] *= 2 * M_PI / LBOX;
	}

	//cout << "maxDistance=" << maxDistance << endl;
	//cout << "nodePointSpacing=" << nodePointSpacing << endl;

	pairDistribution.name = "pairDistribution";
	//pairDistribution.InitGrid(0.0, halfLength, 0.05);
	//pairDistribution.InitGrid(0.0, nodes[nodes.size() - 1], nodes[nodes.size() - 1] / 250.0);
	pairDistribution.InitGrid(0.0, nodes[nodes.size() - 1], nodes[nodes.size() - 1] / (nodes.size() - 2));
	pairDistribution.InitScaling();
	pairDistribution.InitObservables( { "g_2(r_ij)" });

	structureFactor.name = "structureFactor";
	structureFactor.InitGrid(kNorms);
	structureFactor.InitObservables( { "S(k)" });

	additionalObservables.Add(&pairDistribution);
	additionalObservables.Add(&structureFactor);
}

void BosonsDiscrete::CalculateOtherLocalOperators(vector<vector<double> >& R)
{
	vector<double> vecrni(DIM);
	vector<double> evecrni(DIM);
	int bin;
	double rni;

	double potentialExtern = 0;
	double potentialExternComplex = 0;
	double potentialIntern = 0;
	double potentialInternComplex = 0;

	//potential parameters (a -> width; b -> height)
	double a = params[0];
	double b = params[1];
	if (params.size() > 2 && this->time >= params[2])
	{
		a = params[3];
		b = params[4];
	}

	//ClearVector (splineSumsD);
	//ClearVector (splineSumsD2);
	ClearVector(splineSumsFDF);
	ClearVector(splineSumsFDB);
	ClearVector(otherLocalOperators);

	for (int n = 0; n < N; n++)
	{
		//external potential energy
		potentialExtern += GetExternalPotential(R[n]);

		for (int i = 0; i < N; i++)
		{
			rni = VectorDisplacementNIC_DIM(R[n], R[i], vecrni);

			if (rni < maxDistance)
			{
				//internal potential energy
				if (i < n)
				{
					//Gauss
					//double rnia = rni / a;
					//potentialIntern += b * exp(-(rnia * rnia) / 2.0);

					//square well
					if (rni < a)
					{
						potentialIntern += b;
					}

					//Rydberg U_0 / (1 + (r/R_0)^6))
					//R_0 = a, U_0 = b;
					//double rniR0 = rni / a;
					//double rniR03 = rniR0 * rniR0 * rniR0;
					//potentialIntern += b / (1 + rniR03 * rniR03);

					//imaginary part
					//double rniAbsorption = rni - 4.5;
					//double absorptionFactor = 0.0;//-3e-6;
					//if (rniAbsorption > 0)
					//{
					//	potentialInternComplex += absorptionFactor * rniAbsorption * rniAbsorption;
					//}
				}
			}
		}
		//kinetic energy with central difference
		for (int d = 0; d < orderFD; d++)
		{
			splineSumsFDF[n][d] = CalculateWFChange(R, n, 0, d+1);
			splineSumsFDB[n][d] = CalculateWFChange(R, n, 0, -(d+1));
		}
	}
	otherLocalOperators[0] = potentialIntern;
	otherLocalOperators[1] = potentialInternComplex;
	otherLocalOperators[2] = potentialExtern;
	otherLocalOperators[3] = potentialExternComplex;
}

vector<double> BosonsDiscrete::GetCenterOfMass(vector<vector<double> >& R)
{
	vector<double> com = { 0.0, 0.0, 0.0 };
	return com;
}

double BosonsDiscrete::GetExternalPotential(vector<double>& r)
{
	return 0;
}

void BosonsDiscrete::CalculateExpectationValues(vector<vector<double> >& R, vector<double>& uR, vector<double>& uI, double phiR, double phiI)
{
	vector<double> vecrni(DIM);
	int bin;
	double rni;

	double potentialExtern = 0;
	double potentialExternComplex = 0;
	double potentialIntern = 0;
	double potentialInternComplex = 0;

	double kineticR = 0;
	double kineticI = 0;
	double kineticSumR = 0;
	double kineticSumI = 0;
	double kineticSumR_FirstD = 0;
	double kineticSumI_FirstD = 0;

	//potential parameters (a -> width; b -> height)
	double a = params[0];
	double b = params[1];
	if (params.size() > 2 && this->time >= params[2])
	{
		a = params[3];
		b = params[4];
	}

	ClearVector(splineSumsFDF);
	ClearVector(splineSumsFDB);
	ClearVector(otherLocalOperators);

	for (int n = 0; n < N; n++)
	{
		//external potential energy
		potentialExtern += GetExternalPotential(R[n]);

		for (int i = 0; i < N; i++)
		{
			rni = VectorDisplacementNIC_DIM(R[n], R[i], vecrni);

			if (rni < maxDistance)
			{
				//internal potential energy
				if (i < n)
				{
					//square well
					if (rni < a)
					{
						potentialIntern += b;
					}
				}
			}
		}
		//kinetic energy with central difference
		for (int d = 0; d < orderFD; d++)
		{
			splineSumsFDF[n][d] = CalculateWFChange(R, n, 0, d+1) - splineSums;
			splineSumsFDB[n][d] = CalculateWFChange(R, n, 0, -(d+1)) - splineSums;
		}

		auto fdf = splineSumsFDF[n];
		auto fdb = splineSumsFDB[n];
		for (unsigned int f = 0; f < fdf.size(); f++)
		{
			//fdf[f][1] += fdf[f][0];
			//fdf[f][numberOfSplines - 2] += fdf[f][numberOfSplines - 1];
			//fdf[f].erase(fdf[f].begin());
			//fdf[f].pop_back();
			//fdb[f][1] += fdb[f][0];
			//fdb[f][numberOfSplines - 2] += fdb[f][numberOfSplines - 1];
			//fdb[f].erase(fdb[f].begin());
			//fdb[f].pop_back();
		}

		kineticSumR +=  + 1.0 * exp(VectorDotProduct(uR, fdb[0]))
					  	- 2.0
						+ 1.0 * exp(VectorDotProduct(uR, fdf[0]));
		kineticSumI +=  + 1.0 * exp(VectorDotProduct(uI, fdb[0]))
					  	- 2.0
						+ 1.0 * exp(VectorDotProduct(uI, fdf[0]));
		kineticSumR_FirstD +=  - 1.0 / 2.0 * exp(VectorDotProduct(uR, fdb[0]))
					  			+ 1.0 / 2.0 * exp(VectorDotProduct(uR, fdf[0]));
		kineticSumI_FirstD +=  - 1.0 / 2.0 * exp(VectorDotProduct(uI, fdb[0]))
								+ 1.0 / 2.0 * exp(VectorDotProduct(uI, fdf[0]));
		//kineticSumR +=  - 1.0 / 12.0 * exp(VectorDotProduct(uR, fdb[1]))
		//				+ 4.0 / 3.0 * exp(VectorDotProduct(uR, fdb[0]))
		//				- 5.0 / 2.0
		//				+ 4.0 / 3.0 * exp(VectorDotProduct(uR, fdf[0]))
		//				- 1.0 / 12.0 * exp(VectorDotProduct(uR, fdf[1]));
		//kineticSumI +=  - 1.0 / 12.0 * exp(VectorDotProduct(uI, fdb[1]))
		//				+ 4.0 / 3.0 * exp(VectorDotProduct(uI, fdb[0]))
		//				- 5.0 / 2.0
		//				+ 4.0 / 3.0 * exp(VectorDotProduct(uI, fdf[0]))
		//				- 1.0 / 12.0 * exp(VectorDotProduct(uI, fdf[1]));
	}
	otherLocalOperators[0] = potentialIntern;
	otherLocalOperators[1] = potentialInternComplex;
	otherLocalOperators[2] = potentialExtern;
	otherLocalOperators[3] = potentialExternComplex;

	localEnergyR = 0;
	localEnergyI = 0;

	ClearVector(localOperatorsMatrix);
	ClearVector(localOperatorlocalEnergyR);
	ClearVector(localOperatorlocalEnergyI);
	ClearVector(otherExpectationValues);

	double delta = nodes[2] - nodes[1];
	double delta2 = delta * delta;

	kineticR = - 1.0 / delta2 * kineticSumR;
	kineticI = - 1.0 / delta2 * kineticSumI;
	if (IMAGINARY_TIME == 1) //INFO: -1/12 + ... does not sum up to zero if all uI are zero
	{
		kineticI = 0.0;
	}

	localEnergyR = kineticR + otherLocalOperators[0] + otherLocalOperators[2];
	localEnergyI = kineticI + otherLocalOperators[1] + otherLocalOperators[3];

	//cout << "potentialExtern=" << potentialExtern << endl;
	//cout << "potentialIntern=" << potentialIntern << endl;
	//cout << "kineticSumR1=" << kineticSumR1 << endl;
	//cout << "kineticSumI1=" << kineticSumI1 << endl;
	//cout << "kineticSumR2=" << kineticSumR2 << endl;
	//cout << "kineticSumI2=" << kineticSumI2 << endl;
	//cout << "kineticSumR1I1=" << kineticSumR1I1 << endl;
	//cout << "kineticR=" << kineticR << endl;

	//localOperators = splineSums; //INFO: for providing the localOperators to TDVMC.cpp via IPhysicalSystem.GetLocalOperators()

	for (int i = 0; i < N_PARAM; i++)
	{
		//this->localOperators[i] = splineSums[i + 1];
		this->localOperators[i] = splineSums[i];
	}
	//this->localOperators[0] += splineSums[0];
	//this->localOperators[N_PARAM - 1] += splineSums[numberOfSplines - 2];
	//this->localOperators[N_PARAM - 1] += splineSums[numberOfSplines - 1];

	for (int k = 0; k < N_PARAM; k++)
	{
		for (int j = 0; j < N_PARAM; j++)
		{
			localOperatorsMatrix[k][j] = localOperators[k] * localOperators[j];
		}
		localOperatorlocalEnergyR[k] = localOperators[k] * localEnergyR;
		localOperatorlocalEnergyI[k] = localOperators[k] * localEnergyI;
	}

	otherExpectationValues[0] = kineticR;
	otherExpectationValues[1] = otherLocalOperators[0];
	otherExpectationValues[2] = wf;
	otherExpectationValues[3] = exponent;
	otherExpectationValues[4] = kineticSumR_FirstD;
	otherExpectationValues[5] = kineticSumI_FirstD;
	//otherExpectationValues[6] = kineticSumR2;
	//otherExpectationValues[7] = kineticSumI2;
	//otherExpectationValues[8] = kineticSumR1I1;

}

void BosonsDiscrete::CalculateExpectationValues(ICorrelatedSamplingData* sample, vector<double>& uR, vector<double>& uI, double phiR, double phiI)
{
	cout << "not implemented" << endl;
}

void BosonsDiscrete::CalculateAdditionalSystemProperties(vector<vector<double> >& R, vector<double>& uR, vector<double>& uI, double phiR, double phiI)
{
	additionalObservables.ClearValues();
	double rni;
	vector<double> vecrni(DIM);

	//pairDistribution
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < i; j++)
		{
			rni = VectorDisplacementNIC_DIM(R[i], R[j], vecrni);
			if (rni < pairDistribution.grid.max)
			{
				pairDistribution.AddToHistogram(0, rni, 1.0);
			}
		}
	}

	//structureFactor
	vector<double> sumSCos;
	vector<double> sumSSin;
	vector<double> sk;
	InitVector(sumSCos, numOfkValues, 0.0);
	InitVector(sumSSin, numOfkValues, 0.0);
	InitVector(sk, numOfkValues, 0.0);
	for (int i = 0; i < N; i++)
	{
		for (int k = 0; k < numOfkValues; k++)
		{
			for (unsigned int kn = 0; kn < kValues[k].size(); kn++)
			{
				sumSCos[k] += cos(kValues[k][kn][0] * R[i][0] + kValues[k][kn][1] * R[i][1] + kValues[k][kn][2] * R[i][2]);
				sumSSin[k] += sin(kValues[k][kn][0] * R[i][0] + kValues[k][kn][1] * R[i][1] + kValues[k][kn][2] * R[i][2]);
			}
		}
	}
	for (int k = 0; k < numOfkValues; k++)
	{
		sk[k] = (sumSCos[k] * sumSCos[k] + sumSSin[k] * sumSSin[k]) / ((double) (N * kValues[k].size()));
		structureFactor.SetValueAtGridIndex(0, k, sk[k]);
	}
}

void BosonsDiscrete::CalculateWavefunction(vector<vector<double> >& R, vector<double>& uR, vector<double>& uI, double phiR, double phiI)
{
	int bin;
	double rni;
	vector<double> vecrni(DIM);

	ClearVector(splineSums);

	for (int n = 0; n < N; n++)
	{
		for (int i = 0; i < n; i++)
		{
			rni = VectorDisplacementNIC_DIM(R[n], R[i], vecrni);
			if (rni < maxDistance)
			{
				auto interval = lower_bound(nodes.begin(), nodes.end(), rni);
				bin = interval - nodes.begin() - 1;

				splineSums[bin] += 1.0;
			}
			else
			{
				cout << "!!!!" << endl;
			}
		}
	}

	double sum = 0;

	for (int i = 0; i < N_PARAM; i++)
	{
		//sum += uR[i] * splineSums[i + 1];
		sum += uR[i] * splineSums[i];
	}
	//sum += uR[0] * splineSums[0];
	//sum += uR[N_PARAM - 1] * splineSums[numberOfSplines - 2];
	//sum += uR[N_PARAM - 1] * splineSums[numberOfSplines - 1];

	//cout << scientific << setprecision(16) << "sum:" << sum << endl;
	exponent = sum;
	wf = exp(exponent + phiR);
}

void BosonsDiscrete::CalculateWavefunction(ICorrelatedSamplingData* sample, vector<double>& uR, vector<double>& uI, double phiR, double phiI)
{
	cout << "not implemented" << endl;
}

vector<double> BosonsDiscrete::CalculateWFChange(vector<vector<double> >& R, int changedParticleIndex, int coord, int steps)
{
	double sum = 0;
	int bin;
	double rni;
	vector<double> vecrni(DIM);
	vector<double> changedSumsPerBin;

	ClearVector(sumOldPerBin);
	ClearVector(sumNewPerBin);
	InitVector(changedSumsPerBin, sumOldPerBin.size(), 0.0);

	for (int i = 0; i < N; i++)
	{
		if (i != changedParticleIndex)
		{
			rni = VectorDisplacementNIC_DIM(R[i], R[changedParticleIndex], vecrni);
			if (rni < maxDistance)
			{
				auto interval = lower_bound(nodes.begin(), nodes.end(), rni);
				bin = interval - nodes.begin() - 1;
				if (bin >= numberOfSplines)
				{
					cout << "!!!!" << endl;
				}

				sumOldPerBin[bin] += 1.0;
			}
			else
			{
				cout << "!!!" << endl;
			}
			vector<double> changedPos = R[changedParticleIndex];
			//INFO: only 1D with uniform grid
			changedPos += steps * (nodes[2] - nodes[1]);
			rni = VectorDisplacementNIC_DIM(R[i], changedPos, vecrni);
			if (rni < maxDistance)
			{
				auto interval = lower_bound(nodes.begin(), nodes.end(), rni);
				bin = interval - nodes.begin() - 1;

				sumNewPerBin[bin] += 1.0;
			}
			else
			{
				cout << "!!!" << endl;
			}
		}
	}

	for (int i = 0; i < numberOfSplines; i++)
	{
		changedSumsPerBin[i] = max(0.0, splineSums[i] - sumOldPerBin[i] + sumNewPerBin[i]);
	}
	return changedSumsPerBin;
}

void BosonsDiscrete::CalculateWFChange(vector<vector<double> >& R, vector<double>& uR, vector<double>& uI, double phiR, double phiI, int changedParticleIndex, vector<double>& oldPosition)
{
	double sum = 0;
	int bin;
	double rni;
	vector<double> vecrni(DIM);

	ClearVector(sumOldPerBin);
	ClearVector(sumNewPerBin);
	for (int i = 0; i < N; i++)
	{
		if (i != changedParticleIndex)
		{
			rni = VectorDisplacementNIC_DIM(R[i], oldPosition, vecrni);
			if (rni < maxDistance)
			{
				auto interval = lower_bound(nodes.begin(), nodes.end(), rni);
				bin = interval - nodes.begin() - 1;
				if (bin >= numberOfSplines)
				{
					cout << "!!!!" << endl;
				}

				sumOldPerBin[bin] += 1.0;
			}
			else
			{
				cout << "!!!" << endl;
			}
			rni = VectorDisplacementNIC_DIM(R[i], R[changedParticleIndex], vecrni);
			if (rni < maxDistance)
			{
				auto interval = lower_bound(nodes.begin(), nodes.end(), rni);
				bin = interval - nodes.begin() - 1;

				sumNewPerBin[bin] += 1.0;
			}
			else
			{
				cout << "!!!" << endl;
			}
		}
	}

	for (int i = 0; i < numberOfSplines; i++)
	{
		splineSumsNew[i] = max(0.0, splineSums[i] - sumOldPerBin[i] + sumNewPerBin[i]);
	}

	for (int i = 0; i < N_PARAM; i++)
	{
		//sum += uR[i] * splineSumsNew[i + 1];
		sum += uR[i] * splineSumsNew[i];
	}
	//sum += uR[0] * splineSumsNew[0];
	//sum += uR[N_PARAM - 1] * splineSumsNew[numberOfSplines - 2];
	//sum += uR[N_PARAM - 1] * splineSumsNew[numberOfSplines - 1];

	exponentNew = sum;
	wfNew = exp(exponentNew + phiR);
}

double BosonsDiscrete::CalculateWFQuotient(vector<vector<double> >& R, vector<double>& uR, vector<double>& uI, double phiR, double phiI, int changedParticleIndex, vector<double>& oldPosition)
{
	CalculateWFChange(R, uR, uI, phiR, phiI, changedParticleIndex, oldPosition);
	double wfQuotient = exp(2.0 * (exponentNew - exponent));

	//cout << "qu=" << wfQuotient << " (" << exponentNew << " - " << exponent << ")" << endl;

	return wfQuotient;
}

void BosonsDiscrete::AcceptMove()
{
	wf = wfNew;
	exponent = exponentNew;
	splineSums = splineSumsNew;
}

void BosonsDiscrete::InitCorrelatedSamplingData(vector<ICorrelatedSamplingData*>& data)
{
}

void BosonsDiscrete::FillCorrelatedSamplingData(ICorrelatedSamplingData* data)
{
}

}
