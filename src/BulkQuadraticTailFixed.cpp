#include "BulkQuadraticTailFixed.h"

BulkQuadraticTailFixed::BulkQuadraticTailFixed(string configDirectory) : IPhysicalSystem()
{
	this->configDirectory = configDirectory;

	this->USE_NORMALIZATION_AND_PHASE = false;
	this->USE_NIC = true;
	this->USE_MOVE_COM_TO_ZERO = false;
}

void BulkQuadraticTailFixed::InitSystem()
{
	grBinCount = 100;
	numOfkValues = 300;
	grBinStartIndex = 5;
	numOfOtherExpectationValues = grBinStartIndex + grBinCount;
	numOfAdditionalSystemProperties = numOfOtherExpectationValues + numOfkValues;

	rijSplit = 1.9;
	double rijSplit2 = rijSplit * rijSplit;
	double rijSplit3 = rijSplit2 * rijSplit;
	double rijSplit4 = rijSplit3 * rijSplit;

	numberOfSplines = N_PARAM + 1;
	halfLength = LBOX / 2.0;
	nodePointSpacing = rijSplit / (double) (numberOfSplines - 3.0);
	maxDistance = halfLength;
	nodePointSpacing2 = nodePointSpacing * nodePointSpacing;
	double maxDistance2 = maxDistance * maxDistance;

	//cut off
	bcFactors.push_back({
		1.0,
		1.0,
		1.0,
		-1.0
	});
	bcFactors.push_back({
		(2.0 * nodePointSpacing2 + 2.0 * nodePointSpacing * rijSplit + rijSplit2) / rijSplit4,
		(-1.0 * nodePointSpacing2 + rijSplit2) / rijSplit4,
		(2.0 * nodePointSpacing2 - 2.0 * nodePointSpacing * rijSplit + rijSplit2) / rijSplit4,
		-1.0 / maxDistance2
	});
	numberOfSpecialParameters = bcFactors.size();
	numberOfStandardParameters = N_PARAM - numberOfSpecialParameters;

	wf = 0.0;
	exponent = 0.0;
	splineSums.resize(numberOfSplines);
	splineSumsD.resize(numberOfSplines);
	for (auto &k : splineSumsD)
	{
		k.resize(N);
		for (auto &n : k)
		{
			n.resize(DIM);
		}
	}
	splineSumsD2.resize(numberOfSplines);
	for (auto &k : splineSumsD2)
	{
		k.resize(N);
	}
	quadraticSum = 0;
	quadraticSumD.resize(N);
	for (auto &n : quadraticSumD)
	{
		n.resize(DIM);
	}
	ClearVector(quadraticSumD);
	quadraticSumD2.resize(N);
	ClearVector(quadraticSumD2);

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
	exponentNew = 0.0;
	quadraticSumNew = 0.0;
	splineSumsNew.resize(numberOfSplines);
	sumOldPerBin.resize(numberOfSplines);
	sumNewPerBin.resize(numberOfSplines);

	grBins.resize(grBinCount);
	ClearVector(grBins);
	grMaxDistance = halfLength;
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

	kValues = ReadFromFile(this->configDirectory + "kVectors.json", 0);
	for (int k = 0; k < numOfkValues; k++)
	{
		for (unsigned int kn = 0; kn < kValues[k].size(); kn++)
		{
			for (int a = 0; a < DIM; a++)
			{
				kValues[k][kn][a] *= 2 * M_PI / LBOX;
			}
		}
	}

	cout << "maxDistance=" << maxDistance << " nodePointSpacing=" << nodePointSpacing << " LBOX=" << LBOX << endl;
}

vector<double> BulkQuadraticTailFixed::GetCenterOfMass(vector<vector<double> >& R)
{
	vector<double> com = {0.0, 0.0, 0.0};
	return com;
}

double BulkQuadraticTailFixed::GetTailCorrectionKinetic(double b)
{
	double value = 0;
	double maxDistance2 = maxDistance * maxDistance;
	double maxDistance3 = maxDistance2 * maxDistance;
	value = -(4.0 * b * M_PI * (2.0 * b + 3.0 * maxDistance2)) / (3.0 * maxDistance3);
	value *= N;
	return value;
}

double BulkQuadraticTailFixed::GetTailCorrectionPotential()
{
	double value = 0.0;
	//INFO: potential(rijSplit=1.9) is on the order of 10^-77. so no correction is needed
	return value;
}

void BulkQuadraticTailFixed::CalculateExpectationValues(vector<vector<double> >& R, vector<double>& uR, vector<double>& uI, double phiR, double phiI)
{
	vector<double> vecrni(DIM);
	vector<double> evecrni(DIM);
	vector<double> tmp(4);
	double temp;
	double rni;
	double interval;
	int bin;
	double res;
	double res2;

	double potentialExtern = 0;
	double potentialIntern = 0;

	//Gauss potential
	//double a = this->time > 1e-04 ? 0.15 : 0.1;
	//double b = 50.0;
	double a = 0.1;
	double b = 100.0;

	double kineticR = 0;
	double kineticI = 0;
	vector<double> vecKineticSumR1(DIM);
	vector<double> vecKineticSumI1(DIM);
	double kineticSumR1 = 0;
	double kineticSumI1 = 0;
	double kineticSumR1I1 = 0;
	double kineticSumR2 = 0;
	double kineticSumI2 = 0;

	double grBinInterval;
	int grBin;

	ClearVector(splineSumsD);
	ClearVector(splineSumsD2);
	ClearVector(quadraticSumD);
	ClearVector(quadraticSumD2);

	localEnergyR = 0;
	localEnergyI = 0;
	ClearVector(localOperators);
    ClearVector(localOperatorsMatrix);
    ClearVector(localOperatorlocalEnergyR);
    ClearVector(localOperatorlocalEnergyI);
    ClearVector(otherExpectationValues);
    ClearVector(grBins);

	for (int n = 0; n < N; n++)
	{
		//external potential energy
		potentialExtern += 0;

		for (int i = 0; i < N; i++)
		{
			rni = VectorDisplacementNIC(R[n], R[i], vecrni);
			if (rni < maxDistance)
			{
				//internal potential energy
				if (i < n)
				{
					double rnia = rni / a;
					potentialIntern += b * exp(-(rnia * rnia) / 2.0);
				}

				//kinetic energy
				//TODO: maybe it is faster to calculate all localOperators[i][n] before and then loop over this precalculated values
				//		otherwise values are calculated multiple times
				if (i != n) //TODO: also include in (i < n) branch -> then only loop over i < n in "for (int i = 0; i < N; i++)"
				{
					if (rni < rijSplit)
					{
						interval = rni / nodePointSpacing;
						bin = int(interval);
						res = interval - bin;
						res2 = res * res;

						tmp[0] = -1.0 / 2.0 * (1.0 - 2.0 * res + res2);
						tmp[1] = 1.0 / 6.0 * (-12.0 * res + 9.0 * res2);
						tmp[2] = 1.0 / 6.0 * (3.0 + 6.0 * res - 9.0 * res2);
						tmp[3] = 1.0 / 2.0 * res2;
						for (int a = 0; a < DIM; a++)
						{
							evecrni[a] = vecrni[a] / rni;
						}
						for (int a = 0; a < DIM; a++)
						{
							for (int b = 0; b < 4; b++)
							{
								splineSumsD[bin + b][n][a] += tmp[b] * evecrni[a] / nodePointSpacing;
							}
						}
						splineSumsD2[bin][n] += 1.0 / nodePointSpacing2 * (1.0 - res) + 2.0 / (nodePointSpacing * rni) * tmp[0];
						splineSumsD2[bin + 1][n] += 1.0 / nodePointSpacing2 * (1.0 / 6.0 * (-12.0 + 18.0 * res)) + 2.0 / (nodePointSpacing * rni) * tmp[1];
						splineSumsD2[bin + 2][n] += 1.0 / nodePointSpacing2 * (1.0 / 6.0 * (6.0 - 18.0 * res)) + 2.0 / (nodePointSpacing * rni) * tmp[2];
						splineSumsD2[bin + 3][n] += 1.0 / nodePointSpacing2 * (res) + 2.0 / (nodePointSpacing * rni) * tmp[3];
					}
					else
					{
						double rniMinus4 = 1.0 / (rni * rni * rni * rni);
						for (int a = 0; a < DIM; a++)
						{
							quadraticSumD[n][a] += -2.0 * rniMinus4 * vecrni[a];
						}
						quadraticSumD2[n] += 2.0 * rniMinus4;
					}
				}
			}
			//binning for g(r)
			if (i < n && rni < grMaxDistance)
			{
				grBinInterval = rni / grNodePointSpacing;
				grBin = int(grBinInterval);
				//grBins[grBin] += 1.0 / (double)pow(grBin + 1, 2);
				grBins[grBin] += 1.0 / grBinVolumes[grBin];
			}
		}
		for (int a = 0; a < DIM; a++)
		{
			vecKineticSumR1[a] = 0.0;
			vecKineticSumI1[a] = 0.0;
		}

		for (int k = 0; k < numberOfStandardParameters; k++)
		{
			for (int a = 0; a < DIM; a++)
			{
				vecKineticSumR1[a] += uR[k] * splineSumsD[k][n][a];
				vecKineticSumI1[a] += uI[k] * splineSumsD[k][n][a];
			}
			kineticSumR2 += uR[k] * splineSumsD2[k][n];
			kineticSumI2 += uI[k] * splineSumsD2[k][n];
		}
		for (int a = 0; a < DIM; a++)
		{
			temp = (bcFactors[0][0] * splineSumsD[numberOfSplines - 3][n][a] + bcFactors[0][1] * splineSumsD[numberOfSplines - 2][n][a] + bcFactors[0][2] * splineSumsD[numberOfSplines - 1][n][a]);
			vecKineticSumR1[a] += uR[N_PARAM - 2] * temp;
			vecKineticSumI1[a] += uI[N_PARAM - 2] * temp;
			temp = (bcFactors[1][0] * splineSumsD[numberOfSplines - 3][n][a] + bcFactors[1][1] * splineSumsD[numberOfSplines - 2][n][a] + bcFactors[1][2] * splineSumsD[numberOfSplines - 1][n][a] + quadraticSumD[n][a]);
			vecKineticSumR1[a] += uR[N_PARAM - 1] * temp;
			vecKineticSumI1[a] += uI[N_PARAM - 1] * temp;
		}
		temp = (bcFactors[0][0] * splineSumsD2[numberOfSplines - 3][n] + bcFactors[0][1] * splineSumsD2[numberOfSplines - 2][n] + bcFactors[0][2] * splineSumsD2[numberOfSplines - 1][n]);
		kineticSumR2 += uR[N_PARAM - 2] * temp;
		kineticSumI2 += uI[N_PARAM - 2] * temp;
		temp = (bcFactors[1][0] * splineSumsD2[numberOfSplines - 3][n] + bcFactors[1][1] * splineSumsD2[numberOfSplines - 2][n] + bcFactors[1][2] * splineSumsD2[numberOfSplines - 1][n] + quadraticSumD2[n]);
		kineticSumR2 += uR[N_PARAM - 1] * temp;
		kineticSumI2 += uI[N_PARAM - 1] * temp;

		kineticSumR1I1 += 2.0 * VectorDotProduct(vecKineticSumR1, vecKineticSumI1);
		kineticSumR1 += VectorNorm2(vecKineticSumR1);
		kineticSumI1 += VectorNorm2(vecKineticSumI1);
	}

	kineticR = - (kineticSumR1 - kineticSumI1 + kineticSumR2);
	kineticI = - (kineticSumR1I1 + kineticSumI2);

	double kineticTailCorrection = GetTailCorrectionKinetic(uR[N_PARAM - 1]);
	double potentialTailCorrection = GetTailCorrectionPotential();
	//cout << "kinetic TailCorrection=" << kineticTailCorrection << endl;
	//cout << "potential TailCorrection=" << potentialTailCorrection << endl;

	localEnergyR = kineticR + potentialIntern + potentialExtern + kineticTailCorrection + potentialTailCorrection;
	localEnergyI = kineticI;

	//cout << "potentialExtern=" << potentialExtern << endl;
	//cout << "potentialIntern=" << potentialIntern << endl;
	//cout << "kineticSumR1=" << kineticSumR1 << endl;
	//cout << "kineticSumI1=" << kineticSumI1 << endl;
	//cout << "kineticSumR2=" << kineticSumR2 << endl;
	//cout << "kineticSumI2=" << kineticSumI2 << endl;
	//cout << "kineticSumR1I1=" << kineticSumR1I1 << endl;
	//cout << "kineticR=" << kineticR << endl;

	for (int i = 0; i < numberOfStandardParameters; i++)
	{
		localOperators[i] = splineSums[i];
	}
	localOperators[N_PARAM - 2] = (bcFactors[0][0] * splineSums[numberOfSplines - 3] + bcFactors[0][1] * splineSums[numberOfSplines - 2] + bcFactors[0][2] * splineSums[numberOfSplines - 1] + bcFactors[0][3]);
	localOperators[N_PARAM - 1] = (bcFactors[1][0] * splineSums[numberOfSplines - 3] + bcFactors[1][1] * splineSums[numberOfSplines - 2] + bcFactors[1][2] * splineSums[numberOfSplines - 1] + bcFactors[1][3] + quadraticSum);

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
	otherExpectationValues[1] = kineticTailCorrection;
	otherExpectationValues[2] = potentialIntern;
	otherExpectationValues[3] = potentialTailCorrection;
	otherExpectationValues[4] = wf;
	for (int i = 0; i < grBinCount; i++)
	{
		otherExpectationValues[i + grBinStartIndex] = grBins[i];
	}
}

void BulkQuadraticTailFixed::CalculateAdditionalSystemProperties(vector<vector<double> >& R, vector<double>& uR, vector<double>& uI, double phiR, double phiI)
{
    ClearVector(additionalSystemProperties);
	CalculateExpectationValues(R, uR, uI, phiR, phiI);
	for (int i = 0; i < numOfOtherExpectationValues; i++)
	{
		additionalSystemProperties[i] = otherExpectationValues[i];
	}

	vector<double> vecrij(DIM);
	vector<double> sumSCos;
	vector<double> sumSSin;
	vector<double> sk;
	sumSCos.resize(numOfkValues);
	ClearVector(sumSCos);
	sumSSin.resize(numOfkValues);
	ClearVector(sumSSin);
	sk.resize(numOfkValues);
	ClearVector(sk);
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
		sk[k] = (sumSCos[k] * sumSCos[k] + sumSSin[k] * sumSSin[k]) / ((double)(N * kValues[k].size()));
		additionalSystemProperties[numOfOtherExpectationValues + k] = sk[k];
	}
}

void BulkQuadraticTailFixed::CalculateWavefunction(vector<vector<double> >& R, vector<double>& uR, vector<double>& uI, double phiR, double phiI)
{
	double sum = 0;
	double interval;
	int bin;
	double res;
	double res2;
	double res3;
	double rni;
	vector<double> vecrni(DIM);
	double tmp;

	ClearVector(splineSums);
	quadraticSum = 0;

	for (int n = 0; n < N; n++)
	{
		for (int i = 0; i < n; i++)
		{
			rni = VectorDisplacementNIC(R[n], R[i], vecrni);
			if (rni < maxDistance)
			{
				if (rni < rijSplit)
				{
					interval = rni / nodePointSpacing;
					bin = int(interval); //INFO: same as floor(interval) for positive values of "interval" but a little bit faster
					res = interval - bin;
					res2 = res * res;
					res3 = res2 * res;

					tmp = -1.0 / 6.0 * (-1.0 + 3.0 * res - 3.0 * res2 + res3);
					splineSums[bin] += tmp;

					tmp = 1.0 / 6.0 * (4.0 - 6.0 * res2 + 3.0 * res3);;
					splineSums[bin + 1] += tmp;

					tmp = 1.0 / 6.0 * (1.0 + 3.0 * res + 3.0 * res2 - 3.0 * res3);
					splineSums[bin + 2] += tmp;

					tmp = 1.0 / 6.0 * res3;
					splineSums[bin + 3] += tmp;
				}
				else
				{
					quadraticSum += 1.0 / (rni * rni);
				}
			}
		}
	}
	for (int i = 0; i < numberOfStandardParameters; i++)
	{
		sum += uR[i] * splineSums[i];
	}
	sum += uR[N_PARAM - 2] * (bcFactors[0][0] * splineSums[numberOfSplines - 3] + bcFactors[0][1] * splineSums[numberOfSplines - 2] + bcFactors[0][2] * splineSums[numberOfSplines - 1] + bcFactors[0][3]);
	sum += uR[N_PARAM - 1] * (bcFactors[1][0] * splineSums[numberOfSplines - 3] + bcFactors[1][1] * splineSums[numberOfSplines - 2] + bcFactors[1][2] * splineSums[numberOfSplines - 1] + bcFactors[1][3] + quadraticSum);

	wf = exp(sum);
	exponent = sum;
}

void BulkQuadraticTailFixed::CalculateWFChange(vector<vector<double> >& R, vector<double>& uR, vector<double>& uI, double phiR, double phiI, int changedParticleIndex, vector<double>& oldPosition)
{
	double sum = 0;
	double quadraticOld = 0;
	double quadraticNew = 0;
	double interval;
	int bin;
	double res;
	double res2;
	double res3;
	double rni;
	vector<double> vecrni(DIM);

	ClearVector(sumOldPerBin);
	ClearVector(sumNewPerBin);
	for (int i = 0; i < N; i++)
	{
		if (i != changedParticleIndex)
		{
			rni = VectorDisplacementNIC(R[i], oldPosition, vecrni);
			if (rni < maxDistance)
			{
				if (rni < rijSplit)
				{
					interval = rni / nodePointSpacing;
					bin = int(interval);
					res = interval - bin;
					res2 = res * res;
					res3 = res2 * res;

					sumOldPerBin[bin] += -1.0 / 6.0 * (-1.0 + 3.0 * res - 3.0 * res2 + res3);
					sumOldPerBin[bin + 1] += 1.0 / 6.0 * (4.0 - 6.0 * res2 + 3.0 * res3);
					sumOldPerBin[bin + 2] += 1.0 / 6.0 * (1.0 + 3.0 * res + 3.0 * res2 - 3.0 * res3);
					sumOldPerBin[bin + 3] += 1.0 / 6.0 * res3;
				}
				else
				{
					quadraticOld += 1.0 / (rni * rni);
				}
			}
			rni = VectorDisplacementNIC(R[i], R[changedParticleIndex], vecrni);
			if (rni < maxDistance)
			{
				if (rni < rijSplit)
				{
					interval = rni / nodePointSpacing;
					bin = int(interval);
					res = interval - bin;
					res2 = res * res;
					res3 = res2 * res;

					sumNewPerBin[bin] += -1.0 / 6.0 * (-1.0 + 3.0 * res - 3.0 * res2 + res3);
					sumNewPerBin[bin + 1] += 1.0 / 6.0 * (4.0 - 6.0 * res2 + 3.0 * res3);
					sumNewPerBin[bin + 2] += 1.0 / 6.0 * (1.0 + 3.0 * res + 3.0 * res2 - 3.0 * res3);
					sumNewPerBin[bin + 3] += 1.0 / 6.0 * res3;
				}
				else
				{
					quadraticNew += 1.0 / (rni * rni);
				}
			}
		}
	}

	quadraticSumNew = max(0.0, quadraticSum - quadraticOld + quadraticNew);
	for (int i = 0; i < numberOfSplines; i++)
	{
		splineSumsNew[i] = max(0.0, splineSums[i] - sumOldPerBin[i] + sumNewPerBin[i]);
	}

	for (int i = 0; i < numberOfStandardParameters; i++)
	{
		sum += uR[i] * splineSumsNew[i];
	}
	sum += uR[N_PARAM - 2] * (bcFactors[0][0] * splineSumsNew[numberOfSplines - 3] + bcFactors[0][1] * splineSumsNew[numberOfSplines - 2] + bcFactors[0][2] * splineSumsNew[numberOfSplines - 1] + bcFactors[0][3]);
	sum += uR[N_PARAM - 1] * (bcFactors[1][0] * splineSumsNew[numberOfSplines - 3] + bcFactors[1][1] * splineSumsNew[numberOfSplines - 2] + bcFactors[1][2] * splineSumsNew[numberOfSplines - 1] + bcFactors[1][3] + quadraticSumNew);

	wfNew = exp(sum);
	exponentNew = sum;
}

double BulkQuadraticTailFixed::CalculateWFQuotient(vector<vector<double> >& R, vector<double>& uR, vector<double>& uI, double phiR, double phiI, int changedParticleIndex, vector<double>& oldPosition)
{
	CalculateWFChange(R, uR, uI, phiR, phiI, changedParticleIndex, oldPosition);
	double wfQuotient = exp(2.0 * (exponentNew - exponent));

	//cout << "qu=" << wfQuotient << endl;

	return wfQuotient;
}

void BulkQuadraticTailFixed::AcceptMove()
{
	wf = wfNew;
	exponent = exponentNew;
	quadraticSum = quadraticSumNew;
	splineSums = splineSumsNew;
}
