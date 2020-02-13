#include "Bosons1D0th.h"

#include "../SplineFactory.h"

namespace PhysicalSystems
{

Bosons1D0th::Bosons1D0th(vector<double>& params, string configDirectory) :
		IPhysicalSystem(params, configDirectory)
{
	this->USE_NORMALIZATION_AND_PHASE = true;
	this->USE_NIC = true;
	this->USE_MOVE_COM_TO_ZERO = false;

	numberOfSplines = 0;

	halfLength = 0;
	maxDistance = 0;
	numberOfSpecialParametersStart = 0;
	numberOfSpecialParametersEnd = 0;
	numberOfStandardParameters = 0;
	np1 = 0;
	np2 = 0;
	np3 = 0;

	orderFD = 0;

	numOfOtherLocalOperators = 0;
	numOfkValues = 0;

	exponentNew = 0;
	changedParticleIndex = 0;
}

void Bosons1D0th::SetNodes(vector<double> n)
{
	this->nodes = n;
	int periodicNodeCount = 0;
	for (int i = 0; i < periodicNodeCount; i++)
	{
		this->nodes.insert(this->nodes.begin(), -n[i + 1]);
		this->nodes.push_back(2.0 * n[n.size() - 1] - n[n.size() - 2 - i]);
	}
}

void Bosons1D0th::InitSystem()
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
	splineWeights = SplineFactory::GetWeights(nodes, 0);

	//cut off
	SplineFactory::SetBoundaryConditions0_1D_OR_1(nodes, bcFactorsStart);
	SplineFactory::SetBoundaryConditions0_1D_CO_1(nodes, bcFactorsEnd, maxDistance);
	numberOfSpecialParametersStart = bcFactorsStart.size();
	numberOfSpecialParametersEnd = bcFactorsEnd.size();
	numberOfStandardParameters = N_PARAM - numberOfSpecialParametersStart - numberOfSpecialParametersEnd;
	np1 = numberOfSpecialParametersStart;
	np2 = np1 + numberOfStandardParameters;
	np3 = np2 + numberOfSpecialParametersEnd;

	if (N_PARAM != numberOfSplines)
	{
		cout << "!!! WRONG NUMBER OF PARAMETERS !!!" << endl;
	}

	wf = 0.0;
	exponent = 0.0;
	InitVector(splineSums, numberOfSplines, 0.0);
	//InitVector(splineSumsD, numberOfSplines, N, DIM, 0.0);
	//InitVector(splineSumsD2, numberOfSplines, N, 0.0);
	orderFD = 5;
	factorsFD = { -1.0/12.0, 4.0/3.0, -5.0/2.0, 4.0/3.0, -1.0/12.0 };
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
	pairDistribution.InitGrid(0.0, nodes[nodes.size() - 1], nodes[nodes.size() - 1] / (N_PARAM - 1));
	pairDistribution.InitScaling();
	pairDistribution.InitObservables( { "g_2(r_ij)" });

	structureFactor.name = "structureFactor";
	structureFactor.InitGrid(kNorms);
	structureFactor.InitObservables( { "S(k)" });

	additionalObservables.Add(&pairDistribution);
	additionalObservables.Add(&structureFactor);
}

void Bosons1D0th::RefreshLocalOperators()
{
	//INFO: BC
	for (int i = 0; i < np1; i++)
	{
		this->localOperators[i] = (bcFactorsStart[i][0] * splineSums[0] + bcFactorsStart[i][1] * splineSums[1] + bcFactorsStart[i][2] * splineSums[2]);
	}
	int spIdx = 0;
	for (int i = np1; i < np2; i++)
	{
		this->localOperators[i] = splineSums[spIdx];
		spIdx++;
	}
	int idx = 0;
	for (int i = np2; i < np3; i++)
	{
		this->localOperators[i] = (bcFactorsEnd[idx][0] * splineSums[numberOfSplines - 3] + bcFactorsEnd[idx][1] * splineSums[numberOfSplines - 2] + bcFactorsEnd[idx][2] * splineSums[numberOfSplines - 1]);
		idx++;
	}
}

void Bosons1D0th::CalculateLocalOperators(vector<vector<double> >& R)
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

				splineSums[bin] += splineWeights[bin][0][0];
			}
			else
			{
				cout << "!!!!" << endl;
			}
		}
	}
	RefreshLocalOperators();
}

void Bosons1D0th::CalculateOtherLocalOperators(vector<vector<double> >& R)
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

vector<double> Bosons1D0th::GetCenterOfMass(vector<vector<double> >& R)
{
	vector<double> com = { 0.0, 0.0, 0.0 };
	return com;
}

double Bosons1D0th::GetExternalPotential(vector<double>& r)
{
	return 0;
}

void Bosons1D0th::CalculateExpectationValues(vector<double>& O, vector<vector<vector<double> > >& sD, vector<vector<double> >& sD2, vector<double>& otherO, vector<double>& gr, vector<double>& uR, vector<double>& uI, double phiR, double phiI)
{
	double tmp = 0;
	double kineticR = 0;
	double kineticI = 0;
	//vector<double> vecKineticSumR1(DIM);
	//vector<double> vecKineticSumI1(DIM);
	//double kineticSumR1 = 0;
	//double kineticSumI1 = 0;
	//double kineticSumR1I1 = 0;
	//double kineticSumR2 = 0;
	//double kineticSumI2 = 0;
	double kineticSumR = 0;
	double kineticSumI = 0;
	vector<double> tmpSums;

	localEnergyR = 0;
	localEnergyI = 0;

	ClearVector(localOperatorsMatrix);
	ClearVector(localOperatorlocalEnergyR);
	ClearVector(localOperatorlocalEnergyI);
	ClearVector(otherExpectationValues);
	InitVector(tmpSums, 8, 0.0);

	double delta = nodes[2] - nodes[1];
	double delta2 = delta * delta;

	for (int n = 0; n < N; n++)
	{
		ClearVector(tmpSums);
		for (int i = 0; i < N_PARAM; i++)
		{
			tmpSums[0] += uR[i] * (splineSumsFDB[n][1][i] - splineSums[i]);
			tmpSums[1] += uR[i] * (splineSumsFDB[n][0][i] - splineSums[i]);
			tmpSums[2] += uR[i] * (splineSumsFDF[n][0][i] - splineSums[i]);
			tmpSums[3] += uR[i] * (splineSumsFDF[n][1][i] - splineSums[i]);
			tmpSums[4] += uI[i] * (splineSumsFDB[n][1][i] - splineSums[i]);
			tmpSums[5] += uI[i] * (splineSumsFDB[n][0][i] - splineSums[i]);
			tmpSums[6] += uI[i] * (splineSumsFDF[n][0][i] - splineSums[i]);
			tmpSums[7] += uI[i] * (splineSumsFDF[n][1][i] - splineSums[i]);
		}
		kineticSumR += -1.0 / 12.0 * exp(tmpSums[0]) + 4.0 / 3.0 * exp(tmpSums[1]) - 5.0 / 2.0 + 4.0 / 3.0 * exp(tmpSums[2]) - 1.0 / 12.0 * exp(tmpSums[3]);
		kineticSumI += -1.0 / 12.0 * exp(tmpSums[4]) + 4.0 / 3.0 * exp(tmpSums[5]) - 5.0 / 2.0 + 4.0 / 3.0 * exp(tmpSums[6]) - 1.0 / 12.0 * exp(tmpSums[7]);
		//kineticSumR += exp(tmpSums[1]) - 2.0 + exp(tmpSums[2]);
		//kineticSumI += exp(tmpSums[5]) - 2.0 + exp(tmpSums[6]);

		//kineticSumR += 15.0/4.0 - 77.0/6.0 * exp(tmpSumsEnd[0]) + 107.0/6.0 * exp(tmpSumsEnd[1]) - 13.0 * exp(tmpSumsEnd[2]) + 61.0/12.0 * exp(tmpSumsEnd[3]) - 5.0/6.0 * exp(tmpSumsEnd[4]);
		//kineticSumI += 15.0/4.0 - 77.0/6.0 * exp(tmpSumsEnd[5]) + 107.0/6.0 * exp(tmpSumsEnd[6]) - 13.0 * exp(tmpSumsEnd[7]) + 61.0/12.0 * exp(tmpSumsEnd[8]) - 5.0/6.0 * exp(tmpSumsEnd[9]);
	}

	kineticR = - 1.0 / delta2 * kineticSumR;
	kineticI = - 1.0 / delta2 * kineticSumI;
	if (IMAGINARY_TIME == 1) //INFO: -1/12 + ... does not sum up to zero if all uI are zero
	{
		kineticI = 0.0;
	}

	localEnergyR = kineticR + otherO[0] + otherO[2];
	localEnergyI = kineticI + otherO[1];

	//cout << "potentialExtern=" << potentialExtern << endl;
	//cout << "potentialIntern=" << potentialIntern << endl;
	//cout << "kineticSumR1=" << kineticSumR1 << endl;
	//cout << "kineticSumI1=" << kineticSumI1 << endl;
	//cout << "kineticSumR2=" << kineticSumR2 << endl;
	//cout << "kineticSumI2=" << kineticSumI2 << endl;
	//cout << "kineticSumR1I1=" << kineticSumR1I1 << endl;
	//cout << "kineticR=" << kineticR << endl;

	localOperators = O; //INFO: for providing the localOperators to TDVMC.cpp via IPhysicalSystem.GetLocalOperators()
	for (int k = 0; k < N_PARAM; k++)
	{
		for (int j = 0; j < N_PARAM; j++)
		{
			localOperatorsMatrix[k][j] = O[k] * O[j];
		}
		localOperatorlocalEnergyR[k] = O[k] * localEnergyR;
		localOperatorlocalEnergyI[k] = O[k] * localEnergyI;
	}

	otherExpectationValues[0] = kineticR;
	otherExpectationValues[1] = otherO[0];
	otherExpectationValues[2] = wf;
	otherExpectationValues[3] = exponent;
	//otherExpectationValues[4] = kineticSumR1;
	//otherExpectationValues[5] = kineticSumI1;
	//otherExpectationValues[6] = kineticSumR2;
	//otherExpectationValues[7] = kineticSumI2;
	//otherExpectationValues[8] = kineticSumR1I1;
	otherExpectationValues[4] = kineticR;
	otherExpectationValues[5] = kineticI;
}

void Bosons1D0th::CalculateExpectationValues(vector<vector<double> >& R, vector<double>& uR, vector<double>& uI, double phiR, double phiI)
{
	RefreshLocalOperators();
	CalculateOtherLocalOperators(R);
	vector<double> dummy_grBins;
	CalculateExpectationValues(this->localOperators, this->splineSumsD, this->splineSumsD2, this->otherLocalOperators, dummy_grBins, uR, uI, phiR, phiI);
}

void Bosons1D0th::CalculateExpectationValues(ICorrelatedSamplingData* sample, vector<double>& uR, vector<double>& uI, double phiR, double phiI)
{
	CSDataBulkSplines* s = dynamic_cast<CSDataBulkSplines*>(sample);
	CalculateExpectationValues(s->localOperators, s->splineSumsD, s->splineSumsD2, s->otherLocalOperators, s->grBins, uR, uI, phiR, phiI);
}

void Bosons1D0th::CalculateAdditionalSystemProperties(vector<vector<double> >& R, vector<double>& uR, vector<double>& uI, double phiR, double phiI)
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

void Bosons1D0th::CalculateWavefunction(vector<double>& O, vector<double>& uR, vector<double>& uI, double phiR, double phiI)
{
	double sum = 0;

	for (int i = 0; i < N_PARAM; i++)
	{
		sum += uR[i] * O[i];
	}

	//cout << scientific << setprecision(16) << "sum:" << sum << endl;
	exponent = sum;
	wf = exp(exponent + phiR);
}

void Bosons1D0th::CalculateWavefunction(vector<vector<double> >& R, vector<double>& uR, vector<double>& uI, double phiR, double phiI)
{
	CalculateLocalOperators(R);
	CalculateWavefunction(this->localOperators, uR, uI, phiR, phiI);
}

void Bosons1D0th::CalculateWavefunction(ICorrelatedSamplingData* sample, vector<double>& uR, vector<double>& uI, double phiR, double phiI)
{
	CSDataBulkSplines* s = dynamic_cast<CSDataBulkSplines*>(sample);
	CalculateWavefunction(s->localOperators, uR, uI, phiR, phiI);
}

vector<double> Bosons1D0th::CalculateWFChange(vector<vector<double> >& R, int changedParticleIndex, int coord, int steps)
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

				sumOldPerBin[bin] += splineWeights[bin][0][0];
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

				sumNewPerBin[bin] += splineWeights[bin][0][0];
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

void Bosons1D0th::CalculateWFChange(vector<vector<double> >& R, vector<double>& uR, vector<double>& uI, double phiR, double phiI, int changedParticleIndex, vector<double>& oldPosition)
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

				sumOldPerBin[bin] += splineWeights[bin][0][0];
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

				sumNewPerBin[bin] += splineWeights[bin][0][0];
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

	//INFO: BCx
	for (int i = 0; i < np1; i++)
	{
		sum += uR[i] * (bcFactorsStart[i][0] * splineSumsNew[0] + bcFactorsStart[i][1] * splineSumsNew[1] + bcFactorsStart[i][2] * splineSumsNew[2]);
	}
	int spIdx = 0;
	for (int i = np1; i < np2; i++)
	{
		sum += uR[i] * splineSumsNew[spIdx];
		spIdx++;
	}
	int idx = 0;
	for (int i = np2; i < np3; i++)
	{
		sum += uR[i] * (bcFactorsEnd[idx][0] * splineSumsNew[numberOfSplines - 3] + bcFactorsEnd[idx][1] * splineSumsNew[numberOfSplines - 2] + bcFactorsEnd[idx][2] * splineSumsNew[numberOfSplines - 1]);
		idx++;
	}

	exponentNew = sum;
	wfNew = exp(exponentNew + phiR);
}

double Bosons1D0th::CalculateWFQuotient(vector<vector<double> >& R, vector<double>& uR, vector<double>& uI, double phiR, double phiI, int changedParticleIndex, vector<double>& oldPosition)
{
	CalculateWFChange(R, uR, uI, phiR, phiI, changedParticleIndex, oldPosition);
	double wfQuotient = exp(2.0 * (exponentNew - exponent));

	//cout << "qu=" << wfQuotient << " (" << exponentNew << " - " << exponent << ")" << endl;

	return wfQuotient;
}

void Bosons1D0th::AcceptMove()
{
	wf = wfNew;
	exponent = exponentNew;
	splineSums = splineSumsNew;
}

void Bosons1D0th::InitCorrelatedSamplingData(vector<ICorrelatedSamplingData*>& data)
{
}

void Bosons1D0th::FillCorrelatedSamplingData(ICorrelatedSamplingData* data)
{
}

}
