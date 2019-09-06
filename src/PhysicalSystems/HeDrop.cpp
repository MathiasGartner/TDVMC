#include "HeDrop.h"

using namespace std;

namespace PhysicalSystems
{

HeDrop::HeDrop(vector<double>& params, string configDirectory) :
		IPhysicalSystem(params, configDirectory)
{
	this->USE_NORMALIZATION_AND_PHASE = true;
	this->USE_NIC = false;
	this->USE_MOVE_COM_TO_ZERO = true;

	numberOfSplines = 0;
	numberOfShortSplines = 0;
	numberOfLargeSplines = 0;
	mcMillanSum = 0;
	constSum = 0;
	//logSum = 0;
	linearSum = 0;
	rijSplit = 0;
	rijSplineSplit = 0;
	nodePointSpacingShort = 0;
	nodePointSpacingShort2 = 0;
	nodePointSpacingLarge = 0;
	nodePointSpacingLarge2 = 0;
	rijTail = 0;

	grBinCount = 0;
	grBinStartIndex = 0;
	grMaxDistance = 0;
	grNodePointSpacing = 0;
	numOfkValues = 0;
	densityProfileMaxDistance = 0;
	numOfDensityProfileValues = 0;
	densityProfileBinInterval = 0;
	densityProfileBin = 0;
	densityProfileNodePointSpacing = 0;

	mcMillanFactor = 0.0;

	exponentNew = 0;
	mcMillanSumNew = 0;
	constSumNew = 0;
	//logSumNew = 0;
	linearSumNew = 0;

	factorFirstSpline1 = 0;
	factorFirstSpline2 = 0;
	factorSecondSpline1 = 0;
	factorSecondSpline2 = 0;
	factorSecondLastSpline = 0;
	factorSecondLastSplineConst = 0;
	//factorSecondLastSplineLog = 0;
	factorSecondLastSplineLinear = 0;
	factorLastSpline = 0;
	factorLastSplineConst = 0;
	//factorLastSplineLog = 0;
	factorLastSplineLinear = 0;
	factorSecondLastShort = 0;
	factorSecondLastLarge = 0;
	factorLastShort = 0;
	factorLastLarge = 0;
	factorFirstShort = 0;
	factorFirstLarge = 0;
	factorSecondShort = 0;
	factorSecondLarge = 0;
}

void HeDrop::InitSystem()
{
	mcMillanFactor = -4.7;

	grBinCount = 200;
	grBinStartIndex = 3;
	numOfkValues = 84;
	numOfDensityProfileValues = 200;
	numOfOtherExpectationValues = 3 + grBinCount + numOfDensityProfileValues;
	numOfAdditionalSystemProperties = numOfOtherExpectationValues + numOfkValues + 1; //+1: <r^2>

	//rijSplit = 1.75;
	rijSplit = 3.0;
	nodePointSpacingShort = 0.1;
	nodePointSpacingLarge = 0.5;
	//numberOfShortSplines = 83;
	numberOfShortSplines = 70;

	numberOfSplines = N_PARAM + 1 + 2;
	numberOfLargeSplines = numberOfSplines - numberOfShortSplines;
	rijSplineSplit = nodePointSpacingShort * (numberOfShortSplines - 3.0) + rijSplit;
	rijTail = nodePointSpacingLarge * (numberOfLargeSplines - 3.0) + rijSplineSplit;
	nodePointSpacingShort2 = pow(nodePointSpacingShort, 2);
	nodePointSpacingLarge2 = pow(nodePointSpacingLarge, 2);

	//rijSplit = 1.75;
	//rijSplineSplit = 8.6;
	//rijTail = 22.1;

	//numberOfSplines = N_PARAM + 1 + 2;	//+3+3: trailing and leading splines
	//-1: phiR, mcMillan
	//+2: connection between splines with different spacing
	//numberOfShortSplines = 53;
	//numberOfLargeSplines = numberOfSplines - numberOfShortSplines;
	//nodePointSpacingShort = (rijSplineSplit - rijSplit) / (double)(numberOfShortSplines - 3.0);
	//nodePointSpacingShort2 = pow(nodePointSpacingShort, 2);
	//nodePointSpacingLarge = (rijTail - rijSplineSplit) / (double)(numberOfLargeSplines - 3.0);
	//nodePointSpacingLarge2 = pow(nodePointSpacingLarge, 2);

	//mcMillan
	factorFirstSpline1 = -2.0 * mcMillanFactor * nodePointSpacingShort * pow(rijSplit, mcMillanFactor - 1.0);
	factorFirstSpline2 = 1.0;
	factorSecondSpline1 = (mcMillanFactor * nodePointSpacingShort + 3.0 * rijSplit) * pow(rijSplit, mcMillanFactor - 1.0) / 2.0;
	factorSecondSpline2 = -1.0 / 2.0;

	//Exponential
	factorSecondLastSpline = -1.0 / 2.0;
	factorSecondLastSplineConst = 3.0 / 2.0;
	//factorSecondLastSplineLog = 3.0 / 2.0 * log(rijTail) - nodePointSpacingLarge / (2.0 * rijTail);
	factorSecondLastSplineLinear = 3.0 / 2.0 * rijTail - nodePointSpacingLarge / 2.0;
	factorLastSpline = 1.0;
	factorLastSplineConst = 0.0;
	//factorLastSplineLog = 2.0 * nodePointSpacingLarge / rijTail;
	factorLastSplineLinear = 2.0 * nodePointSpacingLarge;

	//Spline connection+
	double d1Plusd2_minus1 = 1.0 / (nodePointSpacingShort + nodePointSpacingLarge);
	factorSecondLastShort = (-nodePointSpacingShort + nodePointSpacingLarge) * d1Plusd2_minus1;
	factorSecondLastLarge = (2.0 * nodePointSpacingLarge) * d1Plusd2_minus1;
	factorLastShort = (-4.0 * nodePointSpacingShort) * d1Plusd2_minus1;
	factorLastLarge = (4.0 * nodePointSpacingLarge) * d1Plusd2_minus1;
	factorFirstShort = (4.0 * nodePointSpacingShort) * d1Plusd2_minus1;
	factorFirstLarge = (-4.0 * nodePointSpacingLarge) * d1Plusd2_minus1;
	factorSecondShort = (2.0 * nodePointSpacingShort) * d1Plusd2_minus1;
	factorSecondLarge = (nodePointSpacingShort - nodePointSpacingLarge) * d1Plusd2_minus1;

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
	constSum = 0;
	constSumD.resize(N);
	for (auto &n : constSumD)
	{
		n.resize(DIM);
	}
	ClearVector(constSumD);
	constSumD2.resize(N);
	ClearVector(constSumD2);
	//logSum = 0;
	//logSumD.resize(N);
	//for (auto &n : logSumD)
	//{
	//	n.resize(DIM);
	//}
	//ClearVector(logSumD);
	//logSumD2.resize(N);
	//ClearVector(logSumD2);
	linearSum = 0;
	linearSumD.resize(N);
	for (auto &n : linearSumD)
	{
		n.resize(DIM);
	}
	ClearVector(linearSumD);
	linearSumD2.resize(N);
	ClearVector(linearSumD2);
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
	constSumNew = 0.0;
	//logSumNew = 0.0;
	linearSumNew = 0.0;
	splineSumsNew.resize(numberOfSplines);
	sumOldPerBin.resize(numberOfSplines);
	sumNewPerBin.resize(numberOfSplines);

	grBins.resize(grBinCount);
	ClearVector(grBins);
	grMaxDistance = rijTail * 2.0;
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
				kValues[k][kn][a] *= 2 * M_PI / LBOX / 2.0;
			}
		}
	}

	densityProfileMaxDistance = grMaxDistance;
	densityProfileNodePointSpacing = grNodePointSpacing;
	densityProfileBins.resize(numOfDensityProfileValues);

	//cout << "nodePointSpacingShort=" << nodePointSpacingShort << ", nodePointSpacingLarge=" << nodePointSpacingLarge << ", rijSplit=" << rijSplit << ", rijSplineSplit=" << rijSplineSplit << ", rijTail=" << rijTail << endl;
}

vector<double> HeDrop::GetCenterOfMass(vector<vector<double> >& R)
{
	vector<double> com = { 0, 0, 0 };
	for (int i = 0; i < N; i++)
	{
		for (int a = 0; a < DIM; a++)
		{
			com[a] += R[i][a];
		}
	}
	for (int a = 0; a < DIM; a++)
	{
		com[a] /= (double) N;
	}
	return com;
}

void HeDrop::CalculateOtherLocalOperators(vector<vector<double> >& R)
{

}

void HeDrop::CalculateExpectationValues(vector<vector<double> >& R, vector<double>& uR, vector<double>& uI, double phiR, double phiI)
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

	double grBinInterval;
	int grBin;

	ClearVector(mcMillanSumD);
	ClearVector(mcMillanSumD2);
	ClearVector(constSumD);
	ClearVector(constSumD2);
	//ClearVector(logSumD);
	//ClearVector(logSumD2);
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

	if (this->time > 1e-6)
	{
		//e = 21.0;
		//rm = 2.35;
	}
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
							//logSumD[n][a] += 1.0 / rni * evecrni[a];
							linearSumD[n][a] += evecrni[a];
						}
						constSumD2[n] += 0.0;
						//logSumD2[n] += pow(rni, -2);
						linearSumD2[n] += 2.0 / rni;
					}
					else
					{
						double nps;
						double nps2;
						if (rni < rijSplineSplit)
						{
							interval = (rni - rijSplit) / nodePointSpacingShort;
							bin = floor(interval);
							res = interval - bin;
							nps = nodePointSpacingShort;
						}
						else
						{
							interval = (rni - rijSplineSplit) / nodePointSpacingLarge;
							bin = floor(interval);
							res = interval - bin;
							bin = bin + numberOfShortSplines;
							nps = nodePointSpacingLarge;
						}
						res2 = pow(res, 2);
						nps2 = pow(nps, 2);

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
								splineSumsD[bin + b][n][a] += tmp[b] * evecrni[a] / nps;
							}
						}
						splineSumsD2[bin][n] += 1.0 / nps2 * (1.0 - res) + 2.0 / (nps * rni) * tmp[0];
						splineSumsD2[bin + 1][n] += 1.0 / nps2 * (1.0 / 6.0 * (-12.0 + 18.0 * res)) + 2.0 / (nps * rni) * tmp[1];
						splineSumsD2[bin + 2][n] += 1.0 / nps2 * (1.0 / 6.0 * (6.0 - 18.0 * res)) + 2.0 / (nps * rni) * tmp[2];
						splineSumsD2[bin + 3][n] += 1.0 / nps2 * (res) + 2.0 / (nps * rni) * tmp[3];
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
		for (int k = 2; k < numberOfShortSplines - 4; k++)
		{
			for (int a = 0; a < DIM; a++)
			{
				vecKineticSumR1[a] += uR[k] * splineSumsD[k + 1][n][a];
				vecKineticSumI1[a] += uI[k] * splineSumsD[k + 1][n][a];
			}
			kineticSumR2 += uR[k] * splineSumsD2[k + 1][n];
			kineticSumI2 += uI[k] * splineSumsD2[k + 1][n];
		}
		temp = (splineSumsD2[numberOfShortSplines - 3][n] + factorSecondLastShort * splineSumsD2[numberOfShortSplines - 1][n] + factorSecondLastLarge * splineSumsD2[numberOfShortSplines][n]);
		kineticSumR2 += uR[numberOfShortSplines - 4] * temp;
		kineticSumI2 += uI[numberOfShortSplines - 4] * temp;
		temp = (splineSumsD2[numberOfShortSplines - 2][n] + factorLastShort * splineSumsD2[numberOfShortSplines - 1][n] + factorLastLarge * splineSumsD2[numberOfShortSplines][n]);
		kineticSumR2 += uR[numberOfShortSplines - 3] * temp;
		kineticSumI2 += uI[numberOfShortSplines - 3] * temp;
		temp = (splineSumsD2[numberOfShortSplines + 1][n] + factorFirstShort * splineSumsD2[numberOfShortSplines - 1][n] + factorFirstLarge * splineSumsD2[numberOfShortSplines][n]);
		kineticSumR2 += uR[numberOfShortSplines - 2] * temp;
		kineticSumI2 += uI[numberOfShortSplines - 2] * temp;
		temp = (splineSumsD2[numberOfShortSplines + 2][n] + factorSecondShort * splineSumsD2[numberOfShortSplines - 1][n] + factorSecondLarge * splineSumsD2[numberOfShortSplines][n]);
		kineticSumR2 += uR[numberOfShortSplines - 1] * temp;
		kineticSumI2 += uI[numberOfShortSplines - 1] * temp;
		for (int a = 0; a < DIM; a++)
		{

			temp = (splineSumsD[numberOfShortSplines - 3][n][a] + factorSecondLastShort * splineSumsD[numberOfShortSplines - 1][n][a] + factorSecondLastLarge * splineSumsD[numberOfShortSplines][n][a]);
			vecKineticSumR1[a] += uR[numberOfShortSplines - 4] * temp;
			vecKineticSumI1[a] += uI[numberOfShortSplines - 4] * temp;
			temp = (splineSumsD[numberOfShortSplines - 2][n][a] + factorLastShort * splineSumsD[numberOfShortSplines - 1][n][a] + factorLastLarge * splineSumsD[numberOfShortSplines][n][a]);
			vecKineticSumR1[a] += uR[numberOfShortSplines - 3] * temp;
			vecKineticSumI1[a] += uI[numberOfShortSplines - 3] * temp;
			temp = (splineSumsD[numberOfShortSplines + 1][n][a] + factorFirstShort * splineSumsD[numberOfShortSplines - 1][n][a] + factorFirstLarge * splineSumsD[numberOfShortSplines][n][a]);
			vecKineticSumR1[a] += uR[numberOfShortSplines - 2] * temp;
			vecKineticSumI1[a] += uI[numberOfShortSplines - 2] * temp;
			temp = (splineSumsD[numberOfShortSplines + 2][n][a] + factorSecondShort * splineSumsD[numberOfShortSplines - 1][n][a] + factorSecondLarge * splineSumsD[numberOfShortSplines][n][a]);
			vecKineticSumR1[a] += uR[numberOfShortSplines - 1] * temp;
			vecKineticSumI1[a] += uI[numberOfShortSplines - 1] * temp;
		}
		for (int k = numberOfShortSplines; k < N_PARAM - 3; k++)
		{
			for (int a = 0; a < DIM; a++)
			{
				vecKineticSumR1[a] += uR[k] * splineSumsD[k + 3][n][a];
				vecKineticSumI1[a] += uI[k] * splineSumsD[k + 3][n][a];
			}
			kineticSumR2 += uR[k] * splineSumsD2[k + 3][n];
			kineticSumI2 += uI[k] * splineSumsD2[k + 3][n];
		}
		for (int a = 0; a < DIM; a++)
		{
			temp = (splineSumsD[numberOfSplines - 3][n][a] + factorSecondLastSpline * splineSumsD[numberOfSplines - 2][n][a] + factorLastSpline * splineSumsD[numberOfSplines - 1][n][a]);
			vecKineticSumR1[a] += uR[N_PARAM - 3] * temp;
			vecKineticSumI1[a] += uI[N_PARAM - 3] * temp;
			temp = (constSumD[n][a] + factorSecondLastSplineConst * splineSumsD[numberOfSplines - 2][n][a] + factorLastSplineConst * splineSumsD[numberOfSplines - 1][n][a]);
			vecKineticSumR1[a] += uR[N_PARAM - 2] * temp;
			vecKineticSumI1[a] += uI[N_PARAM - 2] * temp;
			//temp = (logSumD[n][a] + factorSecondLastSplineLog * splineSumsD[numberOfSplines - 2][n][a] + factorLastSplineLog * splineSumsD[numberOfSplines - 1][n][a]);
			//vecKineticSumR1[a] += uR[N_PARAM - 2] * temp;
			//vecKineticSumI1[a] += uI[N_PARAM - 2] * temp;
			temp = (linearSumD[n][a] + factorSecondLastSplineLinear * splineSumsD[numberOfSplines - 2][n][a] + factorLastSplineLinear * splineSumsD[numberOfSplines - 1][n][a]);
			vecKineticSumR1[a] += uR[N_PARAM - 1] * temp;
			vecKineticSumI1[a] += uI[N_PARAM - 1] * temp;
		}
		temp = (splineSumsD2[numberOfSplines - 3][n] + factorSecondLastSpline * splineSumsD2[numberOfSplines - 2][n] + factorLastSpline * splineSumsD2[numberOfSplines - 1][n]);
		kineticSumR2 += uR[N_PARAM - 3] * temp;
		kineticSumI2 += uI[N_PARAM - 3] * temp;
		temp = (constSumD2[n] + factorSecondLastSplineConst * splineSumsD2[numberOfSplines - 2][n] + factorLastSplineConst * splineSumsD2[numberOfSplines - 1][n]);
		kineticSumR2 += uR[N_PARAM - 2] * temp;
		kineticSumI2 += uI[N_PARAM - 2] * temp;
		//temp = (logSumD2[n] + factorSecondLastSplineLog * splineSumsD2[numberOfSplines - 2][n] + factorLastSplineLog * splineSumsD2[numberOfSplines - 1][n]);
		//kineticSumR2 += uR[N_PARAM - 2] * temp;
		//kineticSumI2 += uI[N_PARAM - 2] * temp;
		temp = (linearSumD2[n] + factorSecondLastSplineLinear * splineSumsD2[numberOfSplines - 2][n] + factorLastSplineLinear * splineSumsD2[numberOfSplines - 1][n]);
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
	for (int i = 2; i < numberOfShortSplines - 4; i++)
	{
		localOperators[i] = splineSums[i + 1];
	}
	localOperators[numberOfShortSplines - 4] = (splineSums[numberOfShortSplines - 3] + factorSecondLastShort * splineSums[numberOfShortSplines - 1] + factorSecondLastLarge * splineSums[numberOfShortSplines]);
	localOperators[numberOfShortSplines - 3] = (splineSums[numberOfShortSplines - 2] + factorLastShort * splineSums[numberOfShortSplines - 1] + factorLastLarge * splineSums[numberOfShortSplines]);
	localOperators[numberOfShortSplines - 2] = (splineSums[numberOfShortSplines + 1] + factorFirstShort * splineSums[numberOfShortSplines - 1] + factorFirstLarge * splineSums[numberOfShortSplines]);
	localOperators[numberOfShortSplines - 1] = (splineSums[numberOfShortSplines + 2] + factorSecondShort * splineSums[numberOfShortSplines - 1] + factorSecondLarge * splineSums[numberOfShortSplines]);
	for (int i = numberOfShortSplines; i < N_PARAM - 3; i++)
	{
		localOperators[i] = splineSums[i + 3];
	}
	localOperators[N_PARAM - 3] = (splineSums[numberOfSplines - 3] + factorSecondLastSpline * splineSums[numberOfSplines - 2] + factorLastSpline * splineSums[numberOfSplines - 1]);
	localOperators[N_PARAM - 2] = (constSum + factorSecondLastSplineConst * splineSums[numberOfSplines - 2] + factorLastSplineConst * splineSums[numberOfSplines - 1]);
	//localOperators[N_PARAM - 2] = (logSum + factorSecondLastSplineLog * splineSums[numberOfSplines - 2] + factorLastSplineLog * splineSums[numberOfSplines - 1]);
	localOperators[N_PARAM - 1] = (linearSum + factorSecondLastSplineLinear * splineSums[numberOfSplines - 2] + factorLastSplineLinear * splineSums[numberOfSplines - 1]);
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

void HeDrop::CalculateExpectationValues(ICorrelatedSamplingData* sample, vector<double>& uR, vector<double>& uI, double phiR, double phiI)
{

}

void HeDrop::CalculateAdditionalSystemProperties(vector<vector<double> >& R, vector<double>& uR, vector<double>& uI, double phiR, double phiI)
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

	double r;
	double r2Sum = 0.0;
	vector<double> com;
	com = GetCenterOfMass(R);
	for (int i = 0; i < N; i++)
	{
		vector<double> diff = R[i] - com;
		r = VectorNorm(diff);
		r2Sum += r * r;
	}
	additionalSystemProperties[numOfAdditionalSystemProperties - 1] = r2Sum / (double) N;
}

void HeDrop::CalculateWavefunction(vector<vector<double> >& R, vector<double>& uR, vector<double>& uI, double phiR, double phiI)
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
	constSum = 0;
	//logSum = 0;
	linearSum = 0;

	for (int n = 0; n < N; n++)
	{
		for (int i = 0; i < n; i++)
		{
			rni = VectorDisplacement(R[n], R[i], vecrni);
			//if (rni < maxDistance)
			{
				if (rni < rijSplit)
				{
					mcMillanSum += pow(rni, mcMillanFactor);
				}
				else if (rni >= rijTail)
				{
					constSum += 1.0;
					//logSum += log(rni);
					linearSum += rni;
				}
				else
				{
					if (rni < rijSplineSplit)
					{
						interval = (rni - rijSplit) / nodePointSpacingShort;
						bin = floor(interval);
						res = interval - bin;
					}
					else
					{
						interval = (rni - rijSplineSplit) / nodePointSpacingLarge;
						bin = floor(interval);
						res = interval - bin;
						bin = bin + numberOfShortSplines;
					}
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
	for (int i = 2; i < numberOfShortSplines - 4; i++)
	{
		sum += uR[i] * splineSums[i + 1];
	}
	sum += uR[numberOfShortSplines - 4] * (splineSums[numberOfShortSplines - 3] + factorSecondLastShort * splineSums[numberOfShortSplines - 1] + factorSecondLastLarge * splineSums[numberOfShortSplines]);
	sum += uR[numberOfShortSplines - 3] * (splineSums[numberOfShortSplines - 2] + factorLastShort * splineSums[numberOfShortSplines - 1] + factorLastLarge * splineSums[numberOfShortSplines]);
	sum += uR[numberOfShortSplines - 2] * (splineSums[numberOfShortSplines + 1] + factorFirstShort * splineSums[numberOfShortSplines - 1] + factorFirstLarge * splineSums[numberOfShortSplines]);
	sum += uR[numberOfShortSplines - 1] * (splineSums[numberOfShortSplines + 2] + factorSecondShort * splineSums[numberOfShortSplines - 1] + factorSecondLarge * splineSums[numberOfShortSplines]);
	for (int i = numberOfShortSplines; i < N_PARAM - 3; i++)
	{
		sum += uR[i] * splineSums[i + 3];
	}
	sum += uR[N_PARAM - 3] * (splineSums[numberOfSplines - 3] + factorSecondLastSpline * splineSums[numberOfSplines - 2] + factorLastSpline * splineSums[numberOfSplines - 1]);
	sum += uR[N_PARAM - 2] * (constSum + factorSecondLastSplineConst * splineSums[numberOfSplines - 2] + factorLastSplineConst * splineSums[numberOfSplines - 1]);
	//sum += uR[N_PARAM - 2] * (logSum + factorSecondLastSplineLog * splineSums[numberOfSplines - 2] + factorLastSplineLog * splineSums[numberOfSplines - 1]);
	sum += uR[N_PARAM - 1] * (linearSum + factorSecondLastSplineLinear * splineSums[numberOfSplines - 2] + factorLastSplineLinear * splineSums[numberOfSplines - 1]);

	wf = exp(sum + phiR);
	exponent = sum;

	//cout << "sum=" << sum << "\t\tsumTmp=" << sumTmp << "\t\twf=" << value << endl;
	//cout << "sum=" << sum << "\t\twf=" << wf << "\t\tlogSum=" << logSum << "\t\tlinearSum=" << linearSum << "\t\tphiR=" << phiR << endl;
}

void HeDrop::CalculateWavefunction(ICorrelatedSamplingData* sample, vector<double>& uR, vector<double>& uI, double phiR, double phiI)
{

}

void HeDrop::CalculateWFChange(vector<vector<double> >& R, vector<double>& uR, vector<double>& uI, double phiR, double phiI, int changedParticleIndex, vector<double>& oldPosition)
{
	double sum = 0;
	double mcMillanOld = 0;
	double mcMillanNew = 0;
	double constOld = 0;
	double constNew = 0;
	//double logOld = 0;
	//double logNew = 0;
	double linearOld = 0;
	double linearNew = 0;
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
			rni = VectorDisplacement(R[i], oldPosition, vecrni);
			//if (rni < maxDistance)
			{
				if (rni < rijSplit)
				{
					mcMillanOld += pow(rni, mcMillanFactor);
				}
				else if (rni >= rijTail)
				{
					constOld += 1.0;
					//logOld += log(rni);
					linearOld += rni;
				}
				else
				{
					if (rni < rijSplineSplit)
					{
						interval = (rni - rijSplit) / nodePointSpacingShort;
						bin = floor(interval);
						res = interval - bin;
					}
					else
					{
						interval = (rni - rijSplineSplit) / nodePointSpacingLarge;
						bin = floor(interval);
						res = interval - bin;
						bin = bin + numberOfShortSplines;
					}
					res2 = pow(res, 2);
					res3 = pow(res, 3);

					sumOldPerBin[bin] += -1.0 / 6.0 * (-1.0 + 3.0 * res - 3.0 * res2 + res3);
					sumOldPerBin[bin + 1] += 1.0 / 6.0 * (4.0 - 6.0 * res2 + 3.0 * res3);
					sumOldPerBin[bin + 2] += 1.0 / 6.0 * (1.0 + 3.0 * res + 3.0 * res2 - 3.0 * res3);
					sumOldPerBin[bin + 3] += 1.0 / 6.0 * res3;
				}
			}
			rni = VectorDisplacement(R[i], R[changedParticleIndex], vecrni);
			//if (rni < maxDistance)
			{
				if (rni < rijSplit)
				{
					mcMillanNew += pow(rni, mcMillanFactor);
				}
				else if (rni >= rijTail)
				{
					constNew += 1.0;
					//logNew += log(rni);
					linearNew += rni;
				}
				else
				{
					if (rni < rijSplineSplit)
					{
						interval = (rni - rijSplit) / nodePointSpacingShort;
						bin = floor(interval);
						res = interval - bin;
					}
					else
					{
						interval = (rni - rijSplineSplit) / nodePointSpacingLarge;
						bin = floor(interval);
						res = interval - bin;
						bin = bin + numberOfShortSplines;
					}
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
	constSumNew = constSum - constOld + constNew;
	//logSumNew = logSum - logOld + logNew;
	linearSumNew = max(0.0, linearSum - linearOld + linearNew);

	sum += uR[0] * (mcMillanSumNew + factorFirstSpline1 * splineSumsNew[0] + factorSecondSpline1 * splineSumsNew[1]);
	sum += uR[1] * (splineSumsNew[2] + factorFirstSpline2 * splineSumsNew[0] + factorSecondSpline2 * splineSumsNew[1]);
	for (int i = 2; i < numberOfShortSplines - 4; i++)
	{
		sum += uR[i] * splineSumsNew[i + 1];
	}
	sum += uR[numberOfShortSplines - 4] * (splineSumsNew[numberOfShortSplines - 3] + factorSecondLastShort * splineSumsNew[numberOfShortSplines - 1] + factorSecondLastLarge * splineSumsNew[numberOfShortSplines]);
	sum += uR[numberOfShortSplines - 3] * (splineSumsNew[numberOfShortSplines - 2] + factorLastShort * splineSumsNew[numberOfShortSplines - 1] + factorLastLarge * splineSumsNew[numberOfShortSplines]);
	sum += uR[numberOfShortSplines - 2] * (splineSumsNew[numberOfShortSplines + 1] + factorFirstShort * splineSumsNew[numberOfShortSplines - 1] + factorFirstLarge * splineSumsNew[numberOfShortSplines]);
	sum += uR[numberOfShortSplines - 1] * (splineSumsNew[numberOfShortSplines + 2] + factorSecondShort * splineSumsNew[numberOfShortSplines - 1] + factorSecondLarge * splineSumsNew[numberOfShortSplines]);
	for (int i = numberOfShortSplines; i < N_PARAM - 3; i++)
	{
		sum += uR[i] * splineSumsNew[i + 3];
	}
	sum += uR[N_PARAM - 3] * (splineSumsNew[numberOfSplines - 3] + factorSecondLastSpline * splineSumsNew[numberOfSplines - 2] + factorLastSpline * splineSumsNew[numberOfSplines - 1]);
	sum += uR[N_PARAM - 2] * (constSumNew + factorSecondLastSplineConst * splineSumsNew[numberOfSplines - 2] + factorLastSplineConst * splineSumsNew[numberOfSplines - 1]);
	//sum += uR[N_PARAM - 2] * (logSumNew + factorSecondLastSplineLog * splineSumsNew[numberOfSplines - 2] + factorLastSplineLog * splineSumsNew[numberOfSplines - 1]);
	sum += uR[N_PARAM - 1] * (linearSumNew + factorSecondLastSplineLinear * splineSumsNew[numberOfSplines - 2] + factorLastSplineLinear * splineSumsNew[numberOfSplines - 1]);

	wfNew = exp(sum + phiR);
	exponentNew = sum;
}

double HeDrop::CalculateWFQuotient(vector<vector<double> >& R, vector<double>& uR, vector<double>& uI, double phiR, double phiI, int changedParticleIndex, vector<double>& oldPosition)
{
	CalculateWFChange(R, uR, uI, phiR, phiI, changedParticleIndex, oldPosition);
	double wfQuotient = exp(2.0 * (exponentNew - exponent));

	//cout << "qu=" << wfQuotient << endl;

	return wfQuotient;
}

void HeDrop::AcceptMove()
{
	wf = wfNew;
	exponent = exponentNew;
	mcMillanSum = mcMillanSumNew;
	for (int i = 0; i < numberOfSplines; i++)
	{
		splineSums[i] = splineSumsNew[i];
	}
	constSum = constSumNew;
	//logSum = logSumNew;
	linearSum = linearSumNew;
}

void HeDrop::InitCorrelatedSamplingData(vector<ICorrelatedSamplingData*>& data)
{
}

void HeDrop::FillCorrelatedSamplingData(ICorrelatedSamplingData* data)
{
}

}

