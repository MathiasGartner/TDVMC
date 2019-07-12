#include "BosonClusterWithLog.h"

#include "../SplineFactory.h"

using namespace std;

namespace PhysicalSystems
{

BosonClusterWithLog::BosonClusterWithLog(vector<double>& params, string configDirectory) :
		BosonCluster(params, configDirectory)
{
	logSum = 0;
	logSumNew = 0;

	logPrefactor = -1.0 / 12.0 / 10.0;
}

void BosonClusterWithLog::InitSystem()
{
	BosonCluster::InitSystem();

	logSum = 0;
	logSumD.resize(N);
	for (auto &n : logSumD)
	{
		n.resize(DIM);
	}
	ClearVector(logSumD);
	logSumD2.resize(N);
	ClearVector(logSumD2);

	logSumNew = 0.0;
}

void BosonClusterWithLog::CalculateExpectationValues(vector<vector<double> >& R, vector<double>& uR, vector<double>& uI, double phiR, double phiI)
{
	double temp;

	int s = 4;
	vector<double> vecrni(DIM);
	vector<double> evecrni(DIM);
	vector<double> tmp1(s);
	vector<double> tmp2(s);
	int bin;
	double rni;
	double rni2;

	double grBinInterval;
	int grBin;

	double potentialExtern = 0;
	double potentialIntern = 0;

	//double sigma = 2.556;
	//double eps = 10.22;
	//double r_V = sigma * 2.5;
	//double r_L = sigma * 2.0;

	//Aziz potential
	//double e = 10.948;
	//double rm = 2.963;
	//double a = 184431.01;
	//double alpha = 10.43329537;
	//double beta = -2.27965105;
	//double d = 1.4826;
	//double c6 = 1.36745214;
	//double c8 = 0.42123807;
	//double c10 = 0.17473318;

	double kineticR = 0;
	double kineticI = 0;
	vector<double> vecKineticSumR1(DIM);
	vector<double> vecKineticSumI1(DIM);
	double kineticSumR1 = 0;
	double kineticSumI1 = 0;
	double kineticSumR1I1 = 0;
	double kineticSumR2 = 0;
	double kineticSumI2 = 0;

	ClearVector(mcMillanSumD);
	ClearVector(mcMillanSumD2);
	ClearVector(constSumD);
	ClearVector(constSumD2);
	ClearVector(logSumD);
	ClearVector(logSumD2);
	ClearVector(linearSumD);
	ClearVector(linearSumD2);
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
	ClearVector(densityProfileBins);
	vector<double> com = GetCenterOfMass(R);

	//if (this->time > 1e-6)
	//{
	//	e = 21.0;
	//	rm = 2.35;
	//}
	for (int n = 0; n < N; n++)
	{
		//external potential energy
		potentialExtern += 0;
		//potentialExtern += pow(sin(10000.0 * this->time), 2) * (pow(R[n][0], 2) + pow(R[n][1], 2) + pow(R[n][2], 2));

		for (int i = 0; i < N; i++)
		{
			rni = VectorDisplacement(R[n], R[i], vecrni);
			//if (rni < maxDistance)
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
					//double rmFactor = 0.0;
					//double eFactor = 0.0;
					//double x = rni / (rm * (1.0 + rmFactor * pow(sin(6000.0 * this->time), 2)));
					//double x2 = pow(x, 2);
					//double xminus2 = 1.0 / x2;
					//double xminus6 = pow(xminus2, 3);
					//double F = 1;
					//if (x < d)
					//{
					//	F = exp(-pow(d / x - 1, 2));
					//}
					//potentialIntern += (1.0 + eFactor * pow(sin(6000.0 * this->time), 2)) * e * (a * exp(-alpha * x + beta * x2) - F * xminus6 * (c6 + xminus2 * (c8 + xminus2 * c10)));

					double sigma = 4.0;
					double eps = 3.56;
					double LJ;
					double sigma_r_6 = pow(sigma / rni, 6);
					LJ = 4.0 * eps * sigma_r_6 * (sigma_r_6 - 1.0);
					potentialIntern += LJ;
				}

				//kinetic energy
				//TODO: maybe it is faster to calculate all localOperators[i][n] before and then loop over this precalculated values
				//		otherwise values are calculated multiple times
				if (i != n) //TODO: also include in (i < n) branch -> then only loop over i < n
				{
					if (rni < rijSplit)
					{
						double rniMinusD = pow(rni, mcMillanFactor - 2.0);
						for (int a = 0; a < DIM; a++)
						{
							mcMillanSumD[n][a] += mcMillanFactor * rniMinusD * vecrni[a];
						}
						mcMillanSumD2[n] += mcMillanFactor * (mcMillanFactor + 1.0) * rniMinusD;
					}
					else if (rni >= rijTail)
					{
						for (int a = 0; a < DIM; a++)
						{
							evecrni[a] = vecrni[a] / rni;
						}
						for (int a = 0; a < DIM; a++)
						{
							constSumD[n][a] += 0.0;
							linearSumD[n][a] += evecrni[a];
						}
						constSumD2[n] += 0.0;
						linearSumD2[n] += 2.0 / rni;
					}
					else
					{
						auto interval = lower_bound(nodes.begin(), nodes.end(), rni);
						bin = interval - nodes.begin() - 1;
						rni2 = rni * rni;

						for (int p = 0; p < s; p++)
						{
							//first derivative
							tmp1[3 - p] = splineWeights[bin - p][p][1] + 2.0 * splineWeights[bin - p][p][2] * rni + 3.0 * splineWeights[bin - p][p][3] * rni2;
							//second derivative
							tmp2[3 - p] = 2.0 * splineWeights[bin - p][p][2] + 6.0 * splineWeights[bin - p][p][3] * rni;
						}

						for (int a = 0; a < DIM; a++)
						{
							evecrni[a] = vecrni[a] / rni;
						}
						for (int a = 0; a < DIM; a++)
						{
							for (int b = 0; b < s; b++)
							{
								splineSumsD[bin - b][n][a] += tmp1[3 - b] * evecrni[a];
							}
						}
						for (int b = 0; b < s; b++)
						{
							double secondDerivativeFactor = DIM - 1.0; //INFO: at least correct for 2D and 3D (and 1D?)
							splineSumsD2[bin - b][n] += tmp2[3 - b] + secondDerivativeFactor / rni * tmp1[3 - b];
						}
					}
					for (int a = 0; a < DIM; a++)
					{
						logSumD[n][a] += 1.0 / rni * evecrni[a];
					}
					logSumD2[n] += pow(rni, -2);
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

		//density profile
		double r;
		vector<double> diff = R[n] - com;
		r = VectorNorm(diff);
		if (r < densityProfileMaxDistance)
		{
			densityProfileBinInterval = r / densityProfileNodePointSpacing;
			densityProfileBin = floor(densityProfileBinInterval);
			densityProfileBins[densityProfileBin] += 1.0 / grBinVolumes[densityProfileBin];
		}

		//kinetic energy
		for (int a = 0; a < DIM; a++)
		{
			vecKineticSumR1[a] = 0.0;
			vecKineticSumI1[a] = 0.0;
		}

		for (int a = 0; a < DIM; a++)
		{
			temp = (mcMillanSumD[n][a] + bcFactors[0][0] * splineSumsD[0][n][a] + bcFactors[0][1] * splineSumsD[1][n][a]);
			vecKineticSumR1[a] += uR[0] * temp;
			vecKineticSumI1[a] += uI[0] * temp;
			temp = (splineSumsD[2][n][a] + bcFactors[1][0] * splineSumsD[0][n][a] + bcFactors[1][1] * splineSumsD[1][n][a]);
			vecKineticSumR1[a] += uR[1] * temp;
			vecKineticSumI1[a] += uI[1] * temp;
		}
		temp = (mcMillanSumD2[n] + bcFactors[0][0] * splineSumsD2[0][n] + bcFactors[0][1] * splineSumsD2[1][n]);
		kineticSumR2 += uR[0] * temp;
		kineticSumI2 += uI[0] * temp;
		temp = (splineSumsD2[2][n] + bcFactors[1][0] * splineSumsD2[0][n] + bcFactors[1][1] * splineSumsD2[1][n]);
		kineticSumR2 += uR[1] * temp;
		kineticSumI2 += uI[1] * temp;

		for (int k = 2; k < N_PARAM - 3; k++)
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
			temp = (splineSumsD[numberOfSplines - 3][n][a] + bcFactors[2][0] * splineSumsD[numberOfSplines - 2][n][a] + bcFactors[2][1] * splineSumsD[numberOfSplines - 1][n][a]);
			vecKineticSumR1[a] += uR[N_PARAM - 3] * temp;
			vecKineticSumI1[a] += uI[N_PARAM - 3] * temp;
			temp = (constSumD[n][a] + bcFactors[3][0] * splineSumsD[numberOfSplines - 2][n][a] + bcFactors[3][1] * splineSumsD[numberOfSplines - 1][n][a]);
			vecKineticSumR1[a] += uR[N_PARAM - 2] * temp;
			vecKineticSumI1[a] += uI[N_PARAM - 2] * temp;
			temp = (linearSumD[n][a] + bcFactors[4][0] * splineSumsD[numberOfSplines - 2][n][a] + bcFactors[4][1] * splineSumsD[numberOfSplines - 1][n][a]);
			vecKineticSumR1[a] += uR[N_PARAM - 1] * temp;
			vecKineticSumI1[a] += uI[N_PARAM - 1] * temp;
		}
		temp = (splineSumsD2[numberOfSplines - 3][n] + bcFactors[2][0] * splineSumsD2[numberOfSplines - 2][n] + bcFactors[2][1] * splineSumsD2[numberOfSplines - 1][n]);
		kineticSumR2 += uR[N_PARAM - 3] * temp;
		kineticSumI2 += uI[N_PARAM - 3] * temp;
		temp = (constSumD2[n] + bcFactors[3][0] * splineSumsD2[numberOfSplines - 2][n] + bcFactors[3][1] * splineSumsD2[numberOfSplines - 1][n]);
		kineticSumR2 += uR[N_PARAM - 2] * temp;
		kineticSumI2 += uI[N_PARAM - 2] * temp;
		temp = (linearSumD2[n] + bcFactors[4][0] * splineSumsD2[numberOfSplines - 2][n] + bcFactors[4][1] * splineSumsD2[numberOfSplines - 1][n]);
		kineticSumR2 += uR[N_PARAM - 1] * temp;
		kineticSumI2 += uI[N_PARAM - 1] * temp;

		for (int a = 0; a < DIM; a++)
		{
			vecKineticSumR1[a] += logPrefactor * logSumD[n][a];
		}
		kineticSumR2 += logPrefactor * logSumD2[n];

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

	localOperators[0] = mcMillanSum + bcFactors[0][0] * splineSums[0] + bcFactors[0][1] * splineSums[1];
	localOperators[1] = splineSums[2] + bcFactors[1][0] * splineSums[0] + bcFactors[1][1] * splineSums[1];
	for (int i = 2; i < N_PARAM - 3; i++)
	{
		localOperators[i] = splineSums[i + 1];
	}
	localOperators[N_PARAM - 3] = splineSums[numberOfSplines - 3] + bcFactors[2][0] * splineSums[numberOfSplines - 2] + bcFactors[2][1] * splineSums[numberOfSplines - 1];
	localOperators[N_PARAM - 2] = constSum + bcFactors[3][0] * splineSums[numberOfSplines - 2] + bcFactors[3][1] * splineSums[numberOfSplines - 1];
	//localOperators[N_PARAM - 2] = (logSum + factorSecondLastSplineLog * splineSums[numberOfSplines - 2] + factorLastSplineLog * splineSums[numberOfSplines - 1]);
	localOperators[N_PARAM - 1] = linearSum + bcFactors[4][0] * splineSums[numberOfSplines - 2] + bcFactors[4][1] * splineSums[numberOfSplines - 1];
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
	for (int i = 0; i < numOfDensityProfileValues; i++)
	{
		otherExpectationValues[i + grBinStartIndex + grBinCount] = densityProfileBins[i];
	}
}

void BosonClusterWithLog::CalculateWavefunction(vector<vector<double> >& R, vector<double>& uR, vector<double>& uI, double phiR, double phiI)
{
	double sum = 0;
	int bin;
	double rni;
	double rni2;
	double rni3;
	vector<double> vecrni(DIM);

	ClearVector(splineSums);
	mcMillanSum = 0;
	constSum = 0;
	logSum = 0;
	linearSum = 0;

	for (int n = 0; n < N; n++)
	{
		for (int i = 0; i < n; i++)
		{
			rni = VectorDisplacement(R[n], R[i], vecrni);
			if (rni < rijSplit)
			{
				mcMillanSum += pow(rni, mcMillanFactor);
			}
			else if (rni >= rijTail)
			{
				constSum += 1.0;
				linearSum += rni;
			}
			else
			{
				auto interval = lower_bound(nodes.begin(), nodes.end(), rni);
				bin = interval - nodes.begin() - 1;
				rni2 = rni * rni;
				rni3 = rni2 * rni;

				for (int p = 0; p < 4; p++)
				{
					splineSums[bin - p] += splineWeights[bin - p][p][0] + splineWeights[bin - p][p][1] * rni + splineWeights[bin - p][p][2] * rni2 + splineWeights[bin - p][p][3] * rni3;
				}
			}
			logSum += log(rni);
		}
	}

	//am
	sum += uR[0] * (mcMillanSum + bcFactors[0][0] * splineSums[0] + bcFactors[0][1] * splineSums[1]);
	//a2
	sum += uR[1] * (splineSums[2] + bcFactors[1][0] * splineSums[0] + bcFactors[1][1] * splineSums[1]);
	for (int i = 2; i < N_PARAM - 3; i++)
	{
		sum += uR[i] * splineSums[i + 1];//TODO: check if *new variables are not used here!
	}
	//a(n-2)
	sum += uR[N_PARAM - 3] * (splineSums[numberOfSplines - 3] + bcFactors[2][0] * splineSums[numberOfSplines - 2] + bcFactors[2][1] * splineSums[numberOfSplines - 1]);
	//ac
	sum += uR[N_PARAM - 2] * (constSum + bcFactors[3][0] * splineSums[numberOfSplines - 2] + bcFactors[3][1] * splineSums[numberOfSplines - 1]);
	//alin
	sum += uR[N_PARAM - 1] * (linearSum + bcFactors[4][0] * splineSums[numberOfSplines - 2] + bcFactors[4][1] * splineSums[numberOfSplines - 1]);

	sum += logPrefactor * logSum;

	wf = exp(sum + phiR);
	exponent = sum;

	//cout << "sum=" << sum << "\t\tsumTmp=" << sumTmp << "\t\twf=" << value << endl;
	//cout << "sum=" << sum << "\t\twf=" << wf << "\t\tlogSum=" << logSum << "\t\tlinearSum=" << linearSum << "\t\tphiR=" << phiR << endl;
}

void BosonClusterWithLog::CalculateWFChange(vector<vector<double> >& R, vector<double>& uR, vector<double>& uI, double phiR, double phiI, int changedParticleIndex, vector<double>& oldPosition)
{
	double sum = 0;

	double mcMillanOld = 0;
	double mcMillanNew = 0;

	int bin;
	double rni;
	double rni2;
	double rni3;
	vector<double> vecrni(DIM);

	double constOld = 0;
	double constNew = 0;
	double logOld = 0;
	double logNew = 0;
	double linearOld = 0;
	double linearNew = 0;

	ClearVector(sumOldPerBin);
	ClearVector(sumNewPerBin);
	for (int i = 0; i < N; i++)
	{
		if (i != changedParticleIndex)
		{
			rni = VectorDisplacement(R[i], oldPosition, vecrni);

			if (rni < rijSplit)
			{
				mcMillanOld += pow(rni, mcMillanFactor);
			}
			else if (rni >= rijTail)
			{
				constOld += 1.0;
				linearOld += rni;
			}
			else
			{
				auto interval = lower_bound(nodes.begin(), nodes.end(), rni);
				bin = interval - nodes.begin() - 1;
				rni2 = rni * rni;
				rni3 = rni2 * rni;

				for (int p = 0; p < 4; p++)
				{
					sumOldPerBin[bin - p] += splineWeights[bin - p][p][0] + splineWeights[bin - p][p][1] * rni + splineWeights[bin - p][p][2] * rni2 + splineWeights[bin - p][p][3] * rni3;
				}
			}
			logOld += log(rni);

			rni = VectorDisplacement(R[i], R[changedParticleIndex], vecrni);
			if (rni < rijSplit)
			{
				mcMillanNew += pow(rni, mcMillanFactor);
			}
			else if (rni >= rijTail)
			{
				constNew += 1.0;
				linearNew += rni;
			}
			else
			{
				auto interval = lower_bound(nodes.begin(), nodes.end(), rni);
				bin = interval - nodes.begin() - 1;
				rni2 = rni * rni;
				rni3 = rni2 * rni;

				for (int p = 0; p < 4; p++)
				{
					sumNewPerBin[bin - p] += splineWeights[bin - p][p][0] + splineWeights[bin - p][p][1] * rni + splineWeights[bin - p][p][2] * rni2 + splineWeights[bin - p][p][3] * rni3;
				}
			}
			logNew += log(rni);
		}
	}

	//TODO: scaling factors? as in NUBosonsBulk?

	mcMillanSumNew = max(0.0, mcMillanSum - mcMillanOld + mcMillanNew);
	for (int i = 0; i < numberOfSplines; i++)
	{
		splineSumsNew[i] = max(0.0, splineSums[i] - sumOldPerBin[i] + sumNewPerBin[i]);
	}
	constSumNew = constSum - constOld + constNew;
	logSumNew = logSum - logOld + logNew;
	linearSumNew = max(0.0, linearSum - linearOld + linearNew);

	//am
	sum += uR[0] * (mcMillanSumNew + bcFactors[0][0] * splineSumsNew[0] + bcFactors[0][1] * splineSumsNew[1]);
	//a2
	sum += uR[1] * (splineSumsNew[2] + bcFactors[1][0] * splineSumsNew[0] + bcFactors[1][1] * splineSumsNew[1]);
	for (int i = 2; i < N_PARAM - 3; i++)
	{
		sum += uR[i] * splineSumsNew[i + 1];
	}
	//a(n-2)
	sum += uR[N_PARAM - 3] * (splineSumsNew[numberOfSplines - 3] + bcFactors[2][0] * splineSumsNew[numberOfSplines - 2] + bcFactors[2][1] * splineSumsNew[numberOfSplines - 1]);
	//ac
	sum += uR[N_PARAM - 2] * (constSumNew + bcFactors[3][0] * splineSumsNew[numberOfSplines - 2] + bcFactors[3][1] * splineSumsNew[numberOfSplines - 1]);
	//alin
	sum += uR[N_PARAM - 1] * (linearSumNew + bcFactors[4][0] * splineSumsNew[numberOfSplines - 2] + bcFactors[4][1] * splineSumsNew[numberOfSplines - 1]);

	sum += logPrefactor * logSumNew;

	wfNew = exp(sum + phiR);
	exponentNew = sum;
}

void BosonClusterWithLog::AcceptMove()
{
	BosonCluster::AcceptMove();

	logSum = logSumNew;
}

}
