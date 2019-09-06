#include "BulkOnlySplines.h"

namespace PhysicalSystems
{

BulkOnlySplines::BulkOnlySplines(vector<double>& params, string configDirectory) :
		IPhysicalSystem(params, configDirectory)
{
	this->USE_NORMALIZATION_AND_PHASE = true;
	this->USE_NIC = true;
	this->USE_MOVE_COM_TO_ZERO = false;

	numberOfSplines = 0;

	halfLength = 0;
	maxDistance = 0;
	nodePointSpacing = 0;
	nodePointSpacing2 = 0;
	numberOfSpecialParameters = 0;
	numberOfStandardParameters = 0;

	grBinCount = 0;
	grBinStartIndex = 0;
	grMaxDistance = 0;
	grNodePointSpacing = 0;
	numOfkValues = 0;

	exponentNew = 0;
	changedParticleIndex = 0;
}

void BulkOnlySplines::InitSystem()
{
	grBinCount = 100;
	grBinStartIndex = 3;
	numOfkValues = 300;
	numOfOtherExpectationValues = 3 + grBinCount;
	numOfAdditionalSystemProperties = numOfOtherExpectationValues + numOfkValues;

	numberOfSplines = N_PARAM + 2;
	halfLength = LBOX / 2.0;
	nodePointSpacing = halfLength / (double) (numberOfSplines - 3.0);
	maxDistance = halfLength;
	nodePointSpacing2 = nodePointSpacing * nodePointSpacing;

	//cut off
	//bcFactors.push_back({
	//	1.0,
	//	-1.0 / 2.0,
	//	1.0,
	//	0.0
	//});
	//bcFactors.push_back({
	//	0.0,
	//	-3.0 / 2.0,
	//	0.0,
	//	1.0
	//});
	bcFactors.push_back( { 1.0, 1.0, 1.0, 0.0 });
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
	splineSumsNew.resize(numberOfSplines);
	sumOldPerBin.resize(numberOfSplines);
	sumNewPerBin.resize(numberOfSplines);

	//this->changedParticleIndex = -1;
	//sumPerBinPerParticle.resize(N - 1);
	//for(unsigned int i = 0; i < sumPerBinPerParticle.size(); i++)
	//{
	//	sumPerBinPerParticle[i].resize(i + 1);
	//	for(unsigned int j = 0; j < sumPerBinPerParticle.size(); j++)
	//	{
	//		sumPerBinPerParticle[i][j].resize(numberOfSplines);
	//	}
	//}
	//sumNewPerBinForChangedParticle.resize(N);
	//for (auto &row : sumNewPerBinForChangedParticle)
	//{
	//	row.resize(numberOfSplines);
	//}

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

vector<double> BulkOnlySplines::GetCenterOfMass(vector<vector<double> >& R)
{
	vector<double> com = { 0.0, 0.0, 0.0 };
	return com;
}

void BulkOnlySplines::CalculateOtherLocalOperators(vector<vector<double> >& R)
{

}

void BulkOnlySplines::CalculateExpectationValues(vector<vector<double> >& R, vector<double>& uR, vector<double>& uI, double phiR, double phiI)
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
	double b = 0.0 * 100.0;

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

	localEnergyR = 0;
	localEnergyI = 0;
	ClearVector(localOperators);
	ClearVector(localOperatorsMatrix);
	ClearVector(localOperatorlocalEnergyR);
	ClearVector(localOperatorlocalEnergyI);
	ClearVector(otherExpectationValues);
	ClearVector(grBins);

	splineSumsOnlyD2 = splineSumsD2;

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

					splineSumsOnlyD2[bin][n] += 1.0 / nodePointSpacing2 * (1.0 - res);
					splineSumsOnlyD2[bin + 1][n] += 1.0 / nodePointSpacing2 * (1.0 / 6.0 * (-12.0 + 18.0 * res));
					splineSumsOnlyD2[bin + 2][n] += 1.0 / nodePointSpacing2 * (1.0 / 6.0 * (6.0 - 18.0 * res));
					splineSumsOnlyD2[bin + 3][n] += 1.0 / nodePointSpacing2 * (res);
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
			vecKineticSumR1[a] += uR[N_PARAM - 1] * temp;
			vecKineticSumI1[a] += uI[N_PARAM - 1] * temp;
			//temp = (bcFactors[1][0] * splineSumsD[numberOfSplines - 3][n][a] + bcFactors[1][1] * splineSumsD[numberOfSplines - 2][n][a] + bcFactors[1][2] * splineSumsD[numberOfSplines - 1][n][a]);
			//vecKineticSumR1[a] += uR[N_PARAM - 1] * temp;
			//vecKineticSumI1[a] += uI[N_PARAM - 1] * temp;
		}
		temp = (bcFactors[0][0] * splineSumsD2[numberOfSplines - 3][n] + bcFactors[0][1] * splineSumsD2[numberOfSplines - 2][n] + bcFactors[0][2] * splineSumsD2[numberOfSplines - 1][n]);
		kineticSumR2 += uR[N_PARAM - 1] * temp;
		kineticSumI2 += uI[N_PARAM - 1] * temp;
		//temp = (bcFactors[1][0] * splineSumsD2[numberOfSplines - 3][n] + bcFactors[1][1] * splineSumsD2[numberOfSplines - 2][n] + bcFactors[1][2] * splineSumsD2[numberOfSplines - 1][n]);
		//kineticSumR2 += uR[N_PARAM - 1] * temp;
		//kineticSumI2 += uI[N_PARAM - 1] * temp;

		kineticSumR1I1 += 2.0 * VectorDotProduct(vecKineticSumR1, vecKineticSumI1);
		kineticSumR1 += VectorNorm2(vecKineticSumR1);
		kineticSumI1 += VectorNorm2(vecKineticSumI1);

		//cout << kineticSumR2 << endl;
	}

	/*
	 vector<vector<vector<double> > > allEvec;
	 vector<vector<double> > allRni;
	 vector<vector<vector<double> > > allSplineSums;
	 allRni.resize(N);
	 allEvec.resize(N);
	 for (int i = 0; i < N; i++)
	 {
	 allRni[i].resize(N);
	 allEvec[i].resize(N);
	 }
	 for (int n = 0; n < N; n++)
	 {
	 for (int i = 0; i < N; i++)
	 {
	 if (i == n)
	 {
	 rni = 0;
	 evecrni = {0, 0, 0};
	 }
	 else
	 {
	 rni = VectorDisplacementNIC(R[n], R[i], vecrni);
	 evecrni = vecrni / rni;
	 }
	 allRni[n][i] = rni;
	 allEvec[n][i] = evecrni;
	 }
	 }
	 allSplineSums.resize(numberOfSplines);
	 for (auto &s : allSplineSums)
	 {
	 s.resize(N);
	 for (auto &p : s)
	 {
	 p.resize(N);
	 }
	 }
	 for (int n = 0; n < N; n++)
	 {
	 for (int i = 0; i < N; i++)
	 {
	 if (allRni[n][i] < maxDistance)
	 {
	 interval = allRni[n][i] / nodePointSpacing;
	 bin = int(interval);
	 res = interval - bin;
	 res2 = res * res;

	 tmp[0] = -1.0 / 2.0 * (1.0 - 2.0 * res + res2);
	 tmp[1] = 1.0 / 6.0 * (-12.0 * res + 9.0 * res2);
	 tmp[2] = 1.0 / 6.0 * (3.0 + 6.0 * res - 9.0 * res2);
	 tmp[3] = 1.0 / 2.0 * res2;

	 for (int b = 0; b < 4; b++)
	 {
	 allSplineSums[bin + b][n][i] = tmp[b] / nodePointSpacing;
	 }
	 }
	 }
	 }
	 double newSum = 0;
	 for (int k = 0; k < N_PARAM; k++)
	 {
	 for (int p = 0; p < N_PARAM; p++)
	 {
	 for (int n = 0; n < N; n++)
	 {
	 for (int i = 0; i < N; i++)
	 {
	 for (int j = 0; j < N; j++)
	 {
	 if (i != n && j != n)
	 {
	 newSum += uR[k] * uR[p] * allSplineSums[k][n][i] * allSplineSums[p][n][j] * VectorDotProduct(allEvec[n][i], allEvec[n][j]);
	 }
	 }
	 }
	 }
	 }
	 }

	 cout << "!!!!!!!!" << newSum << ", " << kineticSumR1 << endl;
	 */

	vector<vector<double> > localKineticEnergiesD1;
	vector<double> localKineticEnergiesD2;

	vector<double> sumD1;
	double sumD2;
	for (int s = 0; s < numberOfStandardParameters; s++)
	{
		sumD1 =
		{	0, 0, 0};
		sumD2 = 0;
		for (int p = 0; p < N; p++)
		{
			for (int a = 0; a < DIM; a++)
			{
				sumD1[a] += uR[s] * splineSumsD[s][p][a];
			}
			//sumD2 += uR[s] * splineSumsOnlyD2[s][p];
			sumD2 += splineSumsOnlyD2[s][p];
		}
		localKineticEnergiesD1.push_back(sumD1);
		localKineticEnergiesD2.push_back(sumD2);
	}
	sumD1 =
	{	0, 0, 0};
	sumD2 = 0;
	for (int p = 0; p < N; p++)
	{
		for (int a = 0; a < DIM; a++)
		{
			sumD1[a] += uR[N_PARAM - 1] * (bcFactors[0][0] * splineSumsD[numberOfSplines - 3][p][a] + bcFactors[0][1] * splineSumsD[numberOfSplines - 2][p][a] + bcFactors[0][2] * splineSumsD[numberOfSplines - 1][p][a]);
		}
		sumD2 += uR[N_PARAM - 1] * (bcFactors[0][0] * splineSumsOnlyD2[numberOfSplines - 3][p] + bcFactors[0][1] * splineSumsOnlyD2[numberOfSplines - 2][p] + bcFactors[0][2] * splineSumsOnlyD2[numberOfSplines - 1][p]);
	}
	localKineticEnergiesD1.push_back(sumD1);
	localKineticEnergiesD2.push_back(sumD2);
	allLocalKineticEnergiesD1.push_back(localKineticEnergiesD1);
	allLocalKineticEnergiesD2.push_back(localKineticEnergiesD2);

	allER1.push_back(kineticSumR1);
	allER2.push_back(kineticSumR2);
	vector<double> d1sums = OuterSum(localKineticEnergiesD1);
	allER1new.push_back(VectorNorm2(d1sums));
	allER2new.push_back(Sum(localKineticEnergiesD2));

	kineticR = -(kineticSumR1 - kineticSumI1 + kineticSumR2);
	kineticI = -(kineticSumR1I1 + kineticSumI2);

	localEnergyR = kineticR + potentialIntern + potentialExtern;
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
	localOperators[N_PARAM - 1] = (bcFactors[0][0] * splineSums[numberOfSplines - 3] + bcFactors[0][1] * splineSums[numberOfSplines - 2] + bcFactors[0][2] * splineSums[numberOfSplines - 1] + bcFactors[0][3]);
	//localOperators[N_PARAM - 1] = (bcFactors[1][0] * splineSums[numberOfSplines - 3] + bcFactors[1][1] * splineSums[numberOfSplines - 2] + bcFactors[1][2] * splineSums[numberOfSplines - 1] +bcFactors[1][3]);
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

	//cout << kineticR << endl;
}

void BulkOnlySplines::CalculateExpectationValues(ICorrelatedSamplingData* sample, vector<double>& uR, vector<double>& uI, double phiR, double phiI)
{

}

void BulkOnlySplines::CalculateAdditionalSystemProperties(vector<vector<double> >& R, vector<double>& uR, vector<double>& uI, double phiR, double phiI)
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

void BulkOnlySplines::CalculateWavefunction(vector<vector<double> >& R, vector<double>& uR, vector<double>& uI, double phiR, double phiI)
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
	//ClearVector(sumPerBinPerParticle);

	for (int n = 0; n < N; n++)
	{
		for (int i = 0; i < n; i++)
		{
			rni = VectorDisplacementNIC(R[n], R[i], vecrni);
			if (rni < maxDistance)
			{
				interval = rni / nodePointSpacing;
				bin = int(interval); //INFO: same as floor(interval) for positive values of "interval" but a little bit faster
				res = interval - bin;
				res2 = res * res;
				res3 = res2 * res;

				tmp = -1.0 / 6.0 * (-1.0 + 3.0 * res - 3.0 * res2 + res3);
				splineSums[bin] += tmp;
				//sumPerBinPerParticle[n - 1][i][bin] = tmp;

				tmp = 1.0 / 6.0 * (4.0 - 6.0 * res2 + 3.0 * res3);
				;
				splineSums[bin + 1] += tmp;
				//sumPerBinPerParticle[n - 1][i][bin + 1] = tmp;

				tmp = 1.0 / 6.0 * (1.0 + 3.0 * res + 3.0 * res2 - 3.0 * res3);
				splineSums[bin + 2] += tmp;
				//sumPerBinPerParticle[n - 1][i][bin + 2] = tmp;

				tmp = 1.0 / 6.0 * res3;
				splineSums[bin + 3] += tmp;
				//sumPerBinPerParticle[n - 1][i][bin + 3] = tmp;
			}
		}
	}
	for (int i = 0; i < numberOfStandardParameters; i++)
	{
		sum += uR[i] * splineSums[i];
	}
	sum += uR[N_PARAM - 1] * (bcFactors[0][0] * splineSums[numberOfSplines - 3] + bcFactors[0][1] * splineSums[numberOfSplines - 2] + bcFactors[0][2] * splineSums[numberOfSplines - 1] + bcFactors[0][3]);
	//sum += uR[N_PARAM - 1] * (bcFactors[1][0] * splineSums[numberOfSplines - 3] + bcFactors[1][1] * splineSums[numberOfSplines - 2] + bcFactors[1][2] * splineSums[numberOfSplines - 1] + bcFactors[1][3]);

	wf = exp(sum + phiR);
	exponent = sum;
}

void BulkOnlySplines::CalculateWavefunction(ICorrelatedSamplingData* sample, vector<double>& uR, vector<double>& uI, double phiR, double phiI)
{

}

void BulkOnlySplines::CalculateWFChange2(vector<vector<double> >& R, vector<double>& uR, vector<double>& uI, double phiR, double phiI, int changedParticleIndex, vector<double>& oldPosition)
{
	//double sum = 0;
	//double interval;
	//int bin;
	//double res;
	//double res2;
	//double res3;
	//double rni;
	//vector<double> vecrni(DIM);
	//double tmp;
	//
	//ClearVector(sumOldPerBin);
	//ClearVector(sumNewPerBin);
	//ClearVector(sumNewPerBinForChangedParticle);
	//this->changedParticleIndex = changedParticleIndex;
	//for (int i = 0; i < N; i++)
	//{
	//	if (i != changedParticleIndex)
	//	{
	//		if (i > changedParticleIndex)
	//		{
	//			sumOldPerBin += sumPerBinPerParticle[i - 1][changedParticleIndex];
	//		}
	//		else
	//		{
	//			sumOldPerBin += sumPerBinPerParticle[changedParticleIndex - 1][i];
	//		}
	//
	//		rni = VectorDisplacementNIC(R[i], R[changedParticleIndex], vecrni);
	//		if (rni < maxDistance)
	//		{
	//			interval = rni / nodePointSpacing;
	//			bin = int(interval);
	//			res = interval - bin;
	//			res2 = res * res;
	//			res3 = res2 * res;
	//
	//			tmp = -1.0 / 6.0 * (-1.0 + 3.0 * res - 3.0 * res2 + res3);
	//			sumNewPerBin[bin] += tmp;
	//			sumNewPerBinForChangedParticle[i][bin] = tmp;
	//
	//			tmp = 1.0 / 6.0 * (4.0 - 6.0 * res2 + 3.0 * res3);
	//			sumNewPerBin[bin + 1] += tmp;
	//			sumNewPerBinForChangedParticle[i][bin + 1] = tmp;
	//
	//			tmp = 1.0 / 6.0 * (1.0 + 3.0 * res + 3.0 * res2 - 3.0 * res3);
	//			sumNewPerBin[bin + 2] += tmp;
	//			sumNewPerBinForChangedParticle[i][bin + 2] = tmp;
	//
	//			tmp = 1.0 / 6.0 * res3;
	//			sumNewPerBin[bin + 3] += tmp;
	//			sumNewPerBinForChangedParticle[i][bin + 3] = tmp;
	//		}
	//	}
	//}
	//
	//for (int i = 0; i < numberOfSplines; i++)
	//{
	//	splineSumsNew[i] = max(0.0, splineSums[i] - sumOldPerBin[i] + sumNewPerBin[i]);
	//}
	//
	//for (int i = 0; i < numberOfStandardParameters; i++)
	//{
	//	sum += uR[i] * splineSumsNew[i];
	//}
	//sum += uR[N_PARAM - 2] * (bcFactors[0][0] * splineSumsNew[numberOfSplines - 3] + bcFactors[0][1] * splineSumsNew[numberOfSplines - 2] + bcFactors[0][2] * splineSumsNew[numberOfSplines - 1] + bcFactors[0][3]);
	//sum += uR[N_PARAM - 1] * (bcFactors[1][0] * splineSumsNew[numberOfSplines - 3] + bcFactors[1][1] * splineSumsNew[numberOfSplines - 2] + bcFactors[1][2] * splineSumsNew[numberOfSplines - 1] + bcFactors[1][3]);
	//
	//wfNew = exp(sum);
	//exponentNew = sum;
}

void BulkOnlySplines::CalculateWFChange(vector<vector<double> >& R, vector<double>& uR, vector<double>& uI, double phiR, double phiI, int changedParticleIndex, vector<double>& oldPosition)
{
	double sum = 0;
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
			rni = VectorDisplacementNIC(R[i], R[changedParticleIndex], vecrni);
			if (rni < maxDistance)
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
		}
	}

	for (int i = 0; i < numberOfSplines; i++)
	{
		splineSumsNew[i] = max(0.0, splineSums[i] - sumOldPerBin[i] + sumNewPerBin[i]);
	}

	for (int i = 0; i < numberOfStandardParameters; i++)
	{
		sum += uR[i] * splineSumsNew[i];
	}
	sum += uR[N_PARAM - 1] * (bcFactors[0][0] * splineSumsNew[numberOfSplines - 3] + bcFactors[0][1] * splineSumsNew[numberOfSplines - 2] + bcFactors[0][2] * splineSumsNew[numberOfSplines - 1] + bcFactors[0][3]);
	//sum += uR[N_PARAM - 1] * (bcFactors[1][0] * splineSumsNew[numberOfSplines - 3] + bcFactors[1][1] * splineSumsNew[numberOfSplines - 2] + bcFactors[1][2] * splineSumsNew[numberOfSplines - 1] + bcFactors[1][3]);

	wfNew = exp(sum + phiR);
	exponentNew = sum;
}

double BulkOnlySplines::CalculateWFQuotient(vector<vector<double> >& R, vector<double>& uR, vector<double>& uI, double phiR, double phiI, int changedParticleIndex, vector<double>& oldPosition)
{
	CalculateWFChange(R, uR, uI, phiR, phiI, changedParticleIndex, oldPosition);
	double wfQuotient = exp(2.0 * (exponentNew - exponent));

	//cout << "qu=" << wfQuotient << endl;

	return wfQuotient;
}

void BulkOnlySplines::AcceptMove()
{
	wf = wfNew;
	exponent = exponentNew;
	splineSums = splineSumsNew;
	//for (int i = 0; i < N - 1; i++)
	//{
	//	if (i > this->changedParticleIndex)
	//	{
	//		sumPerBinPerParticle[i - 1][this->changedParticleIndex] = sumNewPerBinForChangedParticle[i];
	//	}
	//	else if(i < this->changedParticleIndex)
	//	{
	//		sumPerBinPerParticle[this->changedParticleIndex - 1][i] = sumNewPerBinForChangedParticle[i];
	//	}
	//}
}

void BulkOnlySplines::InitCorrelatedSamplingData(vector<ICorrelatedSamplingData*>& data)
{
}

void BulkOnlySplines::FillCorrelatedSamplingData(ICorrelatedSamplingData* data)
{
}

}
