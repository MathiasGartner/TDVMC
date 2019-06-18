#include "HeBulk.h"

namespace PhysicalSystems
{

HeBulk::HeBulk(vector<double>& params, string configDirectory) :
		IPhysicalSystem(params, configDirectory)
{
	this->USE_NORMALIZATION_AND_PHASE = false;
	this->USE_NIC = true;
	this->USE_MOVE_COM_TO_ZERO = false;

	numberOfSplines = 0;
	mcMillanSum = 0;
	rijSplit = 0;
	halfLength = 0;
	maxDistance = 0;
	nodePointSpacing = 0;
	nodePointSpacing2 = 0;

	grBinCount = 0;
	grBinStartIndex = 0;
	grMaxDistance = 0;
	grNodePointSpacing = 0;
	numOfkValues = 0;

	exponentNew = 0;
	mcMillanSumNew = 0;

	factorFirstSpline1 = 0;
	factorFirstSpline2 = 0;
	factorSecondSpline1 = 0;
	factorSecondSpline2 = 0;
	factorLastSpline = 0;
	factorSecondLastSpline = 0;
	factorLastSplinePhi = 0;
	factorSecondLastSplinePhi = 0;
}

void HeBulk::InitSystem()
{
	grBinCount = 100;
	grBinStartIndex = 3;
	numOfkValues = 300;
	numOfOtherExpectationValues = 3 + grBinCount;
	numOfAdditionalSystemProperties = numOfOtherExpectationValues + numOfkValues;

	rijSplit = 1.95;

	numberOfSplines = N_PARAM - 1 + 3 + 3;	//+3+3: trailing and leading splines
											//-1: phiR, mcMillan
	halfLength = LBOX / 2.0;
	nodePointSpacing = (halfLength - rijSplit) / (double) (numberOfSplines - 3.0);
	maxDistance = halfLength;
	//nodePointSpacing = (maxDistance - rijSplit) / (double)(numberOfSplines - 3.0 - 3.0);
	nodePointSpacing2 = pow(nodePointSpacing, 2);

	//mcMillan
	factorFirstSpline1 = 10.0 * nodePointSpacing / pow(rijSplit, 6.0);
	factorFirstSpline2 = 1.0;
	factorSecondSpline1 = (-5.0 * nodePointSpacing + 3.0 * rijSplit) / (2.0 * pow(rijSplit, 6.0));
	factorSecondSpline2 = -1.0 / 2.0;

	//cut off
	factorSecondLastSpline = -1.0 / 2.0;
	factorLastSpline = 1.0;
	factorSecondLastSplinePhi = -3.0 / 2.0;
	factorLastSplinePhi = 0.0;

	wf = 0.0;
	exponent = 0.0;
	mcMillanSum = 0;
	mcMillanSumD.resize(N);
	for (auto &n : mcMillanSumD)
	{
		n.resize(DIM);
	}
	ClearVector(mcMillanSumD);
	mcMillanSumD2.resize(N);
	ClearVector(mcMillanSumD2);
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
	mcMillanSumNew = 0.0;
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

	kValues = ReadKValuesFromJsonFile(this->configDirectory + "kVectors.json");
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

	cout << "maxDistance=" << maxDistance << endl;
	cout << "nodePointSpacing=" << nodePointSpacing << endl;
}

vector<double> HeBulk::GetCenterOfMass(vector<vector<double> >& R)
{
	vector<double> com = { 0.0, 0.0, 0.0 };
	return com;
}

void HeBulk::CalculateOtherLocalOperators(vector<vector<double> >& R)
{

}

void HeBulk::CalculateExpectationValues(vector<vector<double> >& R, vector<double>& uR, vector<double>& uI, double phiR, double phiI)
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

	//double sigma = 2.556;
	//double eps = 10.22;
	//double r_V = sigma * 2.5;
	//double r_L = sigma * 2.0;

	//Aziz potential
	double e = 10.948;
	double rm = 2.963;
	double a = 184431.01;
	double alpha = 10.43329537;
	double beta = -2.27965105;
	double d = 1.4826;
	double c6 = 1.36745214;
	double c8 = 0.42123807;
	double c10 = 0.17473318;

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

	ClearVector(mcMillanSumD);
	ClearVector(mcMillanSumD2);
	ClearVector(splineSumsD);
	ClearVector(splineSumsD2);

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
				//if (i < n && rni < r_V) //TODO: combine kinetic and potential energy if-statements
				//{
				//	double LJ;
				//	double sigma_r_6 = pow(sigma / rni, 6);
				//	LJ = 4.0 * eps * sigma_r_6 * (sigma_r_6 - 1.0);
				//	if (rni > r_L)
				//	{
				//		double S;
				//		S = 1.0 - pow(rni - r_L, 2) * (3.0 * r_V - r_L - 2.0 * rni) / pow(r_V - r_L, 3);
				//		LJ = LJ * S;
				//	}
				//	//LJ = min(LJ, 50.0);
				//	potentialIntern += LJ;
				//}
				if (i < n)
				{
					double x = rni / rm;
					double x2 = pow(x, 2);
					double xminus2 = 1.0 / x2;
					double xminus6 = pow(xminus2, 3);
					double F = 1;
					if (x < d)
					{
						F = exp(-pow(d / x - 1, 2));
					}
					potentialIntern += e * (a * exp(-alpha * x + beta * x2) - F * xminus6 * (c6 + xminus2 * (c8 + xminus2 * c10)));
				}

				//kinetic energy
				//TODO: maybe it is faster to calculate all localOperators[i][n] before and then loop over this precalculated values
				//		otherwise values are calculated multiple times
				if (i != n) //TODO: also include in (i < n) branch -> then only loop over i < n
				{
					if (rni < rijSplit)
					{
						double rniMinus7 = pow(rni, -7);
						for (int a = 0; a < DIM; a++)
						{
							mcMillanSumD[n][a] += -5.0 * rniMinus7 * vecrni[a];
						}
						mcMillanSumD2[n] += 20.0 * rniMinus7;
					}
					else
					{
						interval = (rni - rijSplit) / nodePointSpacing;
						bin = floor(interval);
						res = interval - bin;
						res2 = pow(res, 2);

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
				}
			}
			//binning for g(r)
			if (i < n && rni < grMaxDistance)
			{
				grBinInterval = rni / grNodePointSpacing;
				grBin = floor(grBinInterval);
				//grBins[grBin] += 1.0 / (double)pow(grBin + 1, 2);
				grBins[grBin] += 1.0 / grBinVolumes[grBin];
			}
		}
		for (int a = 0; a < DIM; a++)
		{
			vecKineticSumR1[a] = 0.0;
			vecKineticSumI1[a] = 0.0;
		}

		temp = (mcMillanSumD2[n] + factorFirstSpline1 * splineSumsD2[0][n] + factorSecondSpline1 * splineSumsD2[1][n]);
		kineticSumR2 += uR[0] * temp;
		kineticSumI2 += uI[0] * temp;
		temp = (splineSumsD2[2][n] + factorFirstSpline2 * splineSumsD2[0][n] + factorSecondSpline2 * splineSumsD2[1][n]);
		kineticSumR2 += uR[1] * temp;
		kineticSumI2 += uI[1] * temp;
		for (int a = 0; a < DIM; a++)
		{
			temp = (mcMillanSumD[n][a] + factorFirstSpline1 * splineSumsD[0][n][a] + factorSecondSpline1 * splineSumsD[1][n][a]);
			vecKineticSumR1[a] += uR[0] * temp;
			vecKineticSumI1[a] += uI[0] * temp;
			temp = (splineSumsD[2][n][a] + factorFirstSpline2 * splineSumsD[0][n][a] + factorSecondSpline2 * splineSumsD[1][n][a]);
			vecKineticSumR1[a] += uR[1] * temp;
			vecKineticSumI1[a] += uI[1] * temp;
		}
		for (int k = 2; k < N_PARAM - 2; k++)
		{
			for (int a = 0; a < DIM; a++)
			{
				vecKineticSumR1[a] += uR[k] * splineSumsD[k + 1][n][a];
				vecKineticSumI1[a] += uI[k] * splineSumsD[k + 1][n][a];
			}
			kineticSumR2 += uR[k] * splineSumsD2[k + 1][n];
			kineticSumI2 += uI[k] * splineSumsD2[k + 1][n];
		}
		for (int a = 0; a < DIM; a++)
		{
			temp = (splineSumsD[numberOfSplines - 6][n][a] + factorSecondLastSpline * splineSumsD[numberOfSplines - 5][n][a] + factorLastSpline * splineSumsD[numberOfSplines - 4][n][a]);
			vecKineticSumR1[a] += uR[N_PARAM - 2] * temp;
			vecKineticSumI1[a] += uI[N_PARAM - 2] * temp;
			temp = (1 + factorSecondLastSplinePhi * splineSumsD[numberOfSplines - 5][n][a] + factorLastSplinePhi * splineSumsD[numberOfSplines - 4][n][a]);
			vecKineticSumR1[a] += uR[N_PARAM - 1] * temp;
			vecKineticSumI1[a] += uI[N_PARAM - 1] * temp;
		}
		temp = (splineSumsD2[numberOfSplines - 6][n] + factorSecondLastSpline * splineSumsD2[numberOfSplines - 5][n] + factorLastSpline * splineSumsD2[numberOfSplines - 4][n]);
		kineticSumR2 += uR[N_PARAM - 2] * temp;
		kineticSumI2 += uI[N_PARAM - 2] * temp;
		temp = (factorSecondLastSplinePhi * splineSumsD2[numberOfSplines - 5][n] + factorLastSplinePhi * splineSumsD2[numberOfSplines - 4][n]);
		kineticSumR2 += uR[N_PARAM - 1] * temp;
		kineticSumI2 += uI[N_PARAM - 1] * temp;

		kineticSumR1I1 += 2.0 * VectorDotProduct(vecKineticSumR1, vecKineticSumI1);
		kineticSumR1 += VectorNorm2(vecKineticSumR1);
		kineticSumI1 += VectorNorm2(vecKineticSumI1);
	}

	kineticR = -HBAR2_2M * (kineticSumR1 - kineticSumI1 + kineticSumR2);
	kineticI = -HBAR2_2M * (kineticSumR1I1 + kineticSumI2);

	localEnergyR = kineticR + potentialIntern + potentialExtern;
	localEnergyI = kineticI;

	//cout << "potentialExtern=" << potentialExtern << endl;
	//cout << "potentialIntern=" << potentialIntern << endl;
	//cout << "kineticSumR1=" << kineticSumR1 << endl;
	//cout << "kineticSumI1=" << kineticSumI1 << endl;
	//cout << "kineticSumR2=" << kineticSumR2 << endl;
	//cout << "kineticSumI2=" << kineticSumI2 << endl;
	//cout << "kineticSumR1I1=" << kineticSumR1I1 << endl;

	localOperators[0] = mcMillanSum + factorFirstSpline1 * splineSums[0] + factorSecondSpline1 * splineSums[1];
	localOperators[1] = splineSums[2] + factorFirstSpline2 * splineSums[0] + factorSecondSpline2 * splineSums[1];
	for (int i = 2; i < N_PARAM - 2; i++)
	{
		localOperators[i] = splineSums[i + 1];
	}
	localOperators[N_PARAM - 2] = (splineSums[numberOfSplines - 6] + factorSecondLastSpline * splineSums[numberOfSplines - 5] + factorLastSpline * splineSums[numberOfSplines - 4]);
	localOperators[N_PARAM - 1] = (1.0 + factorSecondLastSplinePhi * splineSums[numberOfSplines - 5] + factorLastSplinePhi * splineSums[numberOfSplines - 4]);
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
	otherExpectationValues[1] = potentialIntern;
	otherExpectationValues[2] = wf;
	for (int i = 0; i < grBinCount; i++)
	{
		otherExpectationValues[i + grBinStartIndex] = grBins[i];
	}
}

void HeBulk::CalculateExpectationValues(ICorrelatedSamplingData* sample, vector<double>& uR, vector<double>& uI, double phiR, double phiI)
{

}

void HeBulk::CalculateAdditionalSystemProperties(vector<vector<double> >& R, vector<double>& uR, vector<double>& uI, double phiR, double phiI)
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
		sk[k] = (sumSCos[k] * sumSCos[k] + sumSSin[k] * sumSSin[k]) / ((double) (N * kValues[k].size()));
		additionalSystemProperties[numOfOtherExpectationValues + k] = sk[k];
	}
}

void HeBulk::CalculateWavefunction(vector<vector<double> >& R, vector<double>& uR, vector<double>& uI, double phiR, double phiI)
{
	double sum = 0;
	double interval;
	int bin;
	double res;
	double res2;
	double res3;
	double rni;
	vector<double> vecrni(DIM);

	ClearVector(splineSums);
	mcMillanSum = 0;

	for (int n = 0; n < N; n++)
	{
		for (int i = 0; i < n; i++)
		{
			rni = VectorDisplacementNIC(R[n], R[i], vecrni);
			if (rni < maxDistance)
			{
				if (rni < rijSplit)
				{
					mcMillanSum += pow(rni, -5.0);
				}
				else
				{
					interval = (rni - rijSplit) / nodePointSpacing;
					bin = floor(interval);
					res = interval - bin;
					res2 = pow(res, 2);
					res3 = pow(res, 3);

					splineSums[bin] += -1.0 / 6.0 * (-1.0 + 3.0 * res - 3.0 * res2 + res3);
					splineSums[bin + 1] += 1.0 / 6.0 * (4.0 - 6.0 * res2 + 3.0 * res3);
					splineSums[bin + 2] += 1.0 / 6.0 * (1.0 + 3.0 * res + 3.0 * res2 - 3.0 * res3);
					splineSums[bin + 3] += 1.0 / 6.0 * res3;
				}
			}
		}
	}

	sum += uR[0] * (mcMillanSum + factorFirstSpline1 * splineSums[0] + factorSecondSpline1 * splineSums[1]);
	sum += uR[1] * (splineSums[2] + factorFirstSpline2 * splineSums[0] + factorSecondSpline2 * splineSums[1]);
	for (int i = 2; i < N_PARAM - 2; i++)
	{
		sum += uR[i] * splineSums[i + 1];
	}
	sum += uR[N_PARAM - 2] * (splineSums[numberOfSplines - 6] + factorSecondLastSpline * splineSums[numberOfSplines - 5] + factorLastSpline * splineSums[numberOfSplines - 4]);
	sum += uR[N_PARAM - 1] * (1.0 + factorSecondLastSplinePhi * splineSums[numberOfSplines - 5] + factorLastSplinePhi * splineSums[numberOfSplines - 4]);

	wf = exp(sum);
	exponent = sum;

	//cout << "sum=" << sum << "\t\tsumTmp=" << sumTmp << "\t\twf=" << value << endl;
}

void HeBulk::CalculateWavefunction(ICorrelatedSamplingData* sample, vector<double>& uR, vector<double>& uI, double phiR, double phiI)
{

}

void HeBulk::CalculateWFChange(vector<vector<double> >& R, vector<double>& uR, vector<double>& uI, double phiR, double phiI, int changedParticleIndex, vector<double>& oldPosition)
{
	double sum = 0;
	double mcMillanOld = 0;
	double mcMillanNew = 0;
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
					mcMillanOld += pow(rni, -5.0);
				}
				else
				{
					interval = (rni - rijSplit) / nodePointSpacing;
					bin = floor(interval);
					res = interval - bin;
					res2 = pow(res, 2);
					res3 = pow(res, 3);

					sumOldPerBin[bin] += -1.0 / 6.0 * (-1.0 + 3.0 * res - 3.0 * res2 + res3);
					sumOldPerBin[bin + 1] += 1.0 / 6.0 * (4.0 - 6.0 * res2 + 3.0 * res3);
					sumOldPerBin[bin + 2] += 1.0 / 6.0 * (1.0 + 3.0 * res + 3.0 * res2 - 3.0 * res3);
					sumOldPerBin[bin + 3] += 1.0 / 6.0 * res3;
				}
			}
			rni = VectorDisplacementNIC(R[i], R[changedParticleIndex], vecrni);
			if (rni < maxDistance)
			{
				if (rni < rijSplit)
				{
					mcMillanNew += pow(rni, -5.0);
				}
				else
				{
					interval = (rni - rijSplit) / nodePointSpacing;
					bin = floor(interval);
					res = interval - bin;
					res2 = pow(res, 2);
					res3 = pow(res, 3);

					sumNewPerBin[bin] += -1.0 / 6.0 * (-1.0 + 3.0 * res - 3.0 * res2 + res3);
					sumNewPerBin[bin + 1] += 1.0 / 6.0 * (4.0 - 6.0 * res2 + 3.0 * res3);
					sumNewPerBin[bin + 2] += 1.0 / 6.0 * (1.0 + 3.0 * res + 3.0 * res2 - 3.0 * res3);
					sumNewPerBin[bin + 3] += 1.0 / 6.0 * res3;
				}
			}
		}
	}

	mcMillanSumNew = max(0.0, mcMillanSum - mcMillanOld + mcMillanNew);
	for (int i = 0; i < numberOfSplines; i++)
	{
		splineSumsNew[i] = max(0.0, splineSums[i] - sumOldPerBin[i] + sumNewPerBin[i]);
	}

	sum += uR[0] * (mcMillanSumNew + factorFirstSpline1 * splineSumsNew[0] + factorSecondSpline1 * splineSumsNew[1]);
	sum += uR[1] * (splineSumsNew[2] + factorFirstSpline2 * splineSumsNew[0] + factorSecondSpline2 * splineSumsNew[1]);
	for (int i = 2; i < N_PARAM - 2; i++)
	{
		sum += uR[i] * splineSumsNew[i + 1];
	}
	sum += uR[N_PARAM - 2] * (splineSumsNew[numberOfSplines - 6] + factorSecondLastSpline * splineSumsNew[numberOfSplines - 5] + factorLastSpline * splineSumsNew[numberOfSplines - 4]);
	sum += uR[N_PARAM - 1] * (1.0 + factorSecondLastSplinePhi * splineSumsNew[numberOfSplines - 5] + factorLastSplinePhi * splineSumsNew[numberOfSplines - 4]);

	//wfNew = exp(sum) * exp(phiR);
	wfNew = exp(sum);
	exponentNew = sum;
}

double HeBulk::CalculateWFQuotient(vector<vector<double> >& R, vector<double>& uR, vector<double>& uI, double phiR, double phiI, int changedParticleIndex, vector<double>& oldPosition)
{
	CalculateWFChange(R, uR, uI, phiR, phiI, changedParticleIndex, oldPosition);
	double wfQuotient = exp(2.0 * (exponentNew - exponent));

	//cout << "qu=" << wfQuotient << endl;

	return wfQuotient;
}

void HeBulk::AcceptMove()
{
	wf = wfNew;
	exponent = exponentNew;
	mcMillanSum = mcMillanSumNew;
	for (int i = 0; i < numberOfSplines; i++)
	{
		splineSums[i] = splineSumsNew[i];
	}
}

void HeBulk::FillCorrelatedSamplingData(ICorrelatedSamplingData* data)
{
}

}
