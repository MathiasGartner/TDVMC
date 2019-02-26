#include "PBosonsBulk.h"
#include "SplineFactory.h"

PBosonsBulk::PBosonsBulk(vector<double>& params, string configDirectory) :
		IPhysicalSystem(params, configDirectory)
{
	this->USE_NORMALIZATION_AND_PHASE = true;
	this->USE_NIC = true;
	this->USE_MOVE_COM_TO_ZERO = false;

	numberOfSplines = 0;

	halfLength = 0;
	maxDistance = 0;
	numberOfSpecialParameters = 0;
	numberOfStandardParameters = 0;

	numOfOtherLocalOperators = 0;
	grBinCount = 0;
	grBinStartIndex = 0;
	grMaxDistance = 0;
	grNodePointSpacing = 0;
	numOfkValues = 0;

	exponentNew = 0;
	changedParticleIndex = 0;
}

void PBosonsBulk::SetNodes(vector<double> n)
{
	this->nodes = n;
}

void PBosonsBulk::InitSystem()
{
	numOfOtherLocalOperators = 4; //INFO: potentialIntern, potentialInternComplex, potentialExtern, potentialExternComplex
	otherLocalOperators.resize(numOfOtherLocalOperators);
	grBinCount = 400;
	grBinStartIndex = 3;
	numOfkValues = 300;
	numOfOtherExpectationValues = 3 + grBinCount;
	numOfAdditionalSystemProperties = numOfOtherExpectationValues + numOfkValues;
	//numOfAdditionalSystemProperties = numOfOtherExpectationValues + numOfkValues + numberOfSplines + numberOfSplines * N * DIM + numberOfSplines * N; //expectationValues + s(k) + splines + splinesD + splinesD2

	//INFO: (N_PARAM + 2) for imaginary time with BC1.1; N_PARAM + 2 - 2 for no BC and real time;
	numberOfSplines = N_PARAM + 2;// - 2;
	halfLength = LBOX / 2.0;
	maxDistance = nodes.empty() ? halfLength : nodes[nodes.size() - 4]; //INFO: for third order splines three nodes should be outside the maxDistance boundary

	if (nodes.empty())
	{
		for (int i = -3; i < N_PARAM + 3;  i++)
		{
			nodes.push_back((i * maxDistance) / ((double)N_PARAM - 1));
		}
	}
	splineWeights = SplineFactory::GetWeights(nodes);

	//cut off
	bcFactors.push_back(SplineFactory::GetBoundaryConditions1_1(nodes));
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

	//TODO: calculate scaling factors according to the node structure
	scalingFactors.resize(numberOfSplines);
	for (int i = 0; i < numberOfSplines; i++)
	{
		scalingFactors[i] = 1.0 / pow((i + 1), DIM - 1); //INFO: dV that is assigned to one spline scales with i^2 for 3D and i for 2D, 1 for 1D
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

	grBins.resize(grBinCount);
	ClearVector(grBins);
	grMaxDistance = maxDistance;
	grNodePointSpacing = grMaxDistance / (double) grBinCount;

	grBinVolumes.resize(grBinCount);
	ClearVector(grBinVolumes);
	for (int i = 0; i < grBinCount; i++)
	{
		if (DIM == 3)
		{
			grBinVolumes[i] = 4.0 * M_PI * pow(grNodePointSpacing * (i + 1), 3.0) / 3.0; //INFO: 3D sphere volume
		}
		else if (DIM == 2)
		{
			grBinVolumes[i] = M_PI * pow(grNodePointSpacing * (i + 1), 2.0); //INFO: 2D "sphere" volume
		}
		else if (DIM == 1)
		{
			grBinVolumes[i] = 2.0 * (grNodePointSpacing * (i + 1)); //INFO: 1D "sphere" volume
		}
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

	//cout << "maxDistance=" << maxDistance << endl;
	//cout << "nodePointSpacing=" << nodePointSpacing << endl;
}

void PBosonsBulk::RefreshLocalOperators()
{
	for (int i = 0; i < numberOfStandardParameters; i++)
	{
		this->localOperators[i] = splineSums[i];
	}
	for (int i = 0; i < numberOfSpecialParameters; i++)
	{
		this->localOperators[N_PARAM - (numberOfSpecialParameters - i)] = (bcFactors[i][0] * splineSums[numberOfSplines - 3] + bcFactors[i][1] * splineSums[numberOfSplines - 2] + bcFactors[i][2] * splineSums[numberOfSplines - 1] + bcFactors[i][3]);
	}
}

void PBosonsBulk::CalculateLocalOperators(vector<vector<double> >& R)
{
	int bin;
	double rni;
	double rni2;
	double rni3;
	vector<double> allDisplacements;

	ClearVector(splineSums);

	for (int n = 0; n < N; n++)
	{
		for (int i = 0; i <= n; i++)
		{
			allDisplacements = AllVectorDisplacements(R[n], R[i], maxDistance);
			for (int d = 0; d < allDisplacements.size(); d++)
			{
				rni = allDisplacements[d];
				auto interval = lower_bound(nodes.begin(), nodes.end(), rni);
				bin = interval - nodes.begin() - 1;
				rni2 = rni * rni;
				rni3 = rni2 * rni;

				for (int p = 0; p < 4; p++)
				{
					splineSums[bin - p] += splineWeights[bin - p][p][0] + splineWeights[bin - p][p][1] * rni + splineWeights[bin - p][p][2] * rni2 + splineWeights[bin - p][p][3] * rni3;
				}
			}
		}
	}
	for (int i = 0; i < numberOfSplines; i++)
	{
		splineSums[i] *= scalingFactors[i];
	}
	RefreshLocalOperators();
}

void PBosonsBulk::CalculateOtherLocalOperators(vector<vector<double> >& R)
{
	int s = 4;
	vector<double> allDisplacements;
	vector<vector<double> > eVecs;
	vector<double> tmp1(s);
	vector<double> tmp2(s);
	int bin;
	double rni;
	double rni2;

	double grBinInterval;
	int grBin;

	double potentialExtern = 0;
	double potentialExternComplex = 0;
	double potentialIntern = 0;
	double potentialInternComplex = 0;

	//Gauss potential
	double a = params[0];
	double b = params[1];
	if (params.size() > 2 && this->time >= params[2])
	{
		a = params[3];
		b = params[4];
	}

	ClearVector(splineSumsD);
	ClearVector(splineSumsD2);
	ClearVector(grBins);
	ClearVector(otherLocalOperators);

	for (int n = 0; n < N; n++)
	{
		//external potential energy
		potentialExtern += GetExternalPotential(R[n]);

		for (int i = 0; i < N; i++)
		{
			allDisplacements = AllVectorDisplacements(R[n], R[i], maxDistance, eVecs);
			for (int d = 0; d < allDisplacements.size(); d++)
			{
				rni = allDisplacements[d];
				//internal potential energy
				if (i < n || rni > halfLength)
				{
					double rnia = rni / a;
					potentialIntern += b * exp(-(rnia * rnia) / 2.0);

					//double rniAbsorption = rni - 4.5;
					//double absorptionFactor = 0.0;//-3e-6;
					//if (rniAbsorption > 0)
					//{
					//	potentialInternComplex += absorptionFactor * rniAbsorption * rniAbsorption;
					//}
				}

				//kinetic energy
				//TODO: maybe it is faster to calculate all localOperators[i][n] before and then loop over this precalculated values
				//		otherwise values are calculated multiple times
				if (i != n || rni > halfLength) //TODO: also include in (i < n) branch -> then only loop over i < n in "for (int i = 0; i < N; i++)"
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
						for (int b = 0; b < s; b++)
						{
							splineSumsD[bin - b][n][a] += tmp1[3 - b] * eVecs[d][a];
						}
					}
					for (int b = 0; b < s; b++)
					{
						double secondDerivativeFactor = DIM - 1.0; //INFO: at least correct for 2D and 3D (and 1D?)
						splineSumsD2[bin - b][n] += tmp2[3 - b] + secondDerivativeFactor/ rni * tmp1[3 - b];
					}
				}
			}
			//TODO: also bin for all PBC distances
			//binning for g(r)
			if (i < n && rni < grMaxDistance)
			{
				grBinInterval = rni / grNodePointSpacing;
				grBin = int(grBinInterval);
				//grBins[grBin] += 1.0 / (double)pow(grBin + 1, 2);
				grBins[grBin] += 1.0 / grBinVolumes[grBin];
			}
		}

		for (int f = 0; f < numberOfSplines; f++)
		{
			splineSumsD[f][n] *= scalingFactors[f];
			splineSumsD2[f][n] *= scalingFactors[f];
		}
	}
	otherLocalOperators[0] = potentialIntern;
	otherLocalOperators[1] = potentialInternComplex;
	otherLocalOperators[2] = potentialExtern;
	otherLocalOperators[3] = potentialExternComplex;
}

vector<double> PBosonsBulk::GetCenterOfMass(vector<vector<double> >& R)
{
	vector<double> com = { 0.0, 0.0, 0.0 };
	return com;
}

double PBosonsBulk::GetExternalPotential(vector<double>& r)
{
	return 0;
}

void PBosonsBulk::CalculateExpectationValues(vector<double>& O, vector<vector<vector<double> > >& sD, vector<vector<double> >& sD2, vector<double>& otherO, vector<double>& gr, vector<double>& uR, vector<double>& uI, double phiR, double phiI)
{
	double temp;

	double kineticR = 0;
	double kineticI = 0;
	vector<double> vecKineticSumR1(DIM);
	vector<double> vecKineticSumI1(DIM);
	double kineticSumR1 = 0;
	double kineticSumI1 = 0;
	double kineticSumR1I1 = 0;
	double kineticSumR2 = 0;
	double kineticSumI2 = 0;

	localEnergyR = 0;
	localEnergyI = 0;

	ClearVector(localOperatorsMatrix);
	ClearVector(localOperatorlocalEnergyR);
	ClearVector(localOperatorlocalEnergyI);
	ClearVector(otherExpectationValues);

	for (int n = 0; n < N; n++)
	{
		for (int a = 0; a < DIM; a++)
		{
			vecKineticSumR1[a] = 0.0;
			vecKineticSumI1[a] = 0.0;
		}

		for (int k = 0; k < numberOfStandardParameters; k++)
		{
			for (int a = 0; a < DIM; a++)
			{
				vecKineticSumR1[a] += uR[k] * sD[k][n][a];
				vecKineticSumI1[a] += uI[k] * sD[k][n][a];
			}
			kineticSumR2 += uR[k] * sD2[k][n];
			kineticSumI2 += uI[k] * sD2[k][n];
		}
		for (int i = 0; i < numberOfSpecialParameters; i++)
		{
			for (int a = 0; a < DIM; a++)
			{
				temp = (bcFactors[i][0] * sD[numberOfSplines - 3][n][a] + bcFactors[i][1] * sD[numberOfSplines - 2][n][a] + bcFactors[i][2] * sD[numberOfSplines - 1][n][a]);
				vecKineticSumR1[a] += uR[N_PARAM - (numberOfSpecialParameters - i)] * temp;
				vecKineticSumI1[a] += uI[N_PARAM - 1] * temp;
			}
			temp = (bcFactors[i][0] * sD2[numberOfSplines - 3][n] + bcFactors[i][1] * sD2[numberOfSplines - 2][n] + bcFactors[i][2] * sD2[numberOfSplines - 1][n]);
			kineticSumR2 += uR[N_PARAM - (numberOfSpecialParameters - i)] * temp;
			kineticSumI2 += uI[N_PARAM - 1] * temp;
		}

		kineticSumR1I1 += 2.0 * VectorDotProduct(vecKineticSumR1, vecKineticSumI1);
		kineticSumR1 += VectorNorm2(vecKineticSumR1);
		kineticSumI1 += VectorNorm2(vecKineticSumI1);
	}

	kineticR = -(kineticSumR1 - kineticSumI1 + kineticSumR2);
	kineticI = -(kineticSumR1I1 + kineticSumI2);

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
	for (int i = 0; i < grBinCount; i++)
	{
		otherExpectationValues[i + grBinStartIndex] = gr[i];
	}
}

void PBosonsBulk::CalculateExpectationValues(vector<vector<double> >& R, vector<double>& uR, vector<double>& uI, double phiR, double phiI)
{
	RefreshLocalOperators();
	CalculateOtherLocalOperators(R);
	CalculateExpectationValues(this->localOperators, this->splineSumsD, this->splineSumsD2, this->otherLocalOperators, this->grBins, uR, uI, phiR, phiI);
}

void PBosonsBulk::CalculateExpectationValues(ICorrelatedSamplingData* sample, vector<double>& uR, vector<double>& uI, double phiR, double phiI)
{
	CSDataBulkSplines* s = dynamic_cast<CSDataBulkSplines*>(sample);
	CalculateExpectationValues(s->localOperators, s->splineSumsD, s->splineSumsD2, s->otherLocalOperators, s->grBins, uR, uI, phiR, phiI);
}

void PBosonsBulk::CalculateAdditionalSystemProperties(vector<vector<double> >& R, vector<double>& uR, vector<double>& uI, double phiR, double phiI)
{
	//int index = 0;

	ClearVector(additionalSystemProperties);
	CalculateExpectationValues(R, uR, uI, phiR, phiI);
	for (int i = 0; i < numOfOtherExpectationValues; i++)
	{
		additionalSystemProperties[i] = otherExpectationValues[i];
		//index++;
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
		//index++;
	}

	//for (int i = 0; i < numberOfSplines; i++)
	//{
	//	additionalSystemProperties[index] = splineSums[i];
	//	index++;
	//}
	//
	//for (int k = 0; k < numberOfSplines; k++)
	//{
	//	for (int n = 0; n < N; n++)
	//	{
	//		for (int a = 0; a < DIM; a++)
	//		{
	//			additionalSystemProperties[index] = splineSumsD[k][n][a];
	//			index++;
	//		}
	//	}
	//}
	//
	//for (int k = 0; k < numberOfSplines; k++)
	//{
	//	for (int n = 0; n < N; n++)
	//	{
	//		additionalSystemProperties[index] = splineSumsD2[k][n];
	//		index++;
	//	}
	//}
}

void PBosonsBulk::CalculateWavefunction(vector<double>& O, vector<double>& uR, vector<double>& uI, double phiR, double phiI)
{
	double sum = 0;

	for (int i = 0; i < N_PARAM; i++)
	{
		sum += uR[i] * O[i];
	}

	exponent = sum;
	wf = exp(exponent + phiR);
}

void PBosonsBulk::CalculateWavefunction(vector<vector<double> >& R, vector<double>& uR, vector<double>& uI, double phiR, double phiI)
{
	CalculateLocalOperators(R);
	CalculateWavefunction(this->localOperators, uR, uI, phiR, phiI);
}

void PBosonsBulk::CalculateWavefunction(ICorrelatedSamplingData* sample, vector<double>& uR, vector<double>& uI, double phiR, double phiI)
{
	CSDataBulkSplines* s = dynamic_cast<CSDataBulkSplines*>(sample);
	CalculateWavefunction(s->localOperators, uR, uI, phiR, phiI);
}

void PBosonsBulk::CalculateWFChange(vector<vector<double> >& R, vector<double>& uR, vector<double>& uI, double phiR, double phiI, int changedParticleIndex, vector<double>& oldPosition)
{
	double sum = 0;
	int bin;
	double rni;
	double rni2;
	double rni3;
	vector<double> allDisplacements;

	ClearVector(sumOldPerBin);
	ClearVector(sumNewPerBin);
	for (int i = 0; i < N; i++)
	{
		if (i != changedParticleIndex)
		{
			allDisplacements = AllVectorDisplacements(R[i], oldPosition, maxDistance);
			for (int d = 0; d < allDisplacements.size(); d++)
			{
				rni = allDisplacements[d];
				auto interval = lower_bound(nodes.begin(), nodes.end(), rni);
				bin = interval - nodes.begin() - 1;
				rni2 = rni * rni;
				rni3 = rni2 * rni;

				for (int p = 0; p < 4; p++)
				{
					sumOldPerBin[bin - p] += splineWeights[bin - p][p][0] + splineWeights[bin - p][p][1] * rni + splineWeights[bin - p][p][2] * rni2 + splineWeights[bin - p][p][3] * rni3;
				}
			}
			allDisplacements = AllVectorDisplacements(R[i], R[changedParticleIndex], maxDistance);
			for (int d = 0; d < allDisplacements.size(); d++)
			{
				rni = allDisplacements[d];
				auto interval = lower_bound(nodes.begin(), nodes.end(), rni);
				bin = interval - nodes.begin() - 1;
				rni2 = rni * rni;
				rni3 = rni2 * rni;

				for (int p = 0; p < 4; p++)
				{
					sumNewPerBin[bin - p] += splineWeights[bin - p][p][0] + splineWeights[bin - p][p][1] * rni + splineWeights[bin - p][p][2] * rni2 + splineWeights[bin - p][p][3] * rni3;
				}
			}
		}
	}

	for (int i = 0; i < numberOfSplines; i++)
	{
		sumOldPerBin[i] *= scalingFactors[i];
		sumNewPerBin[i] *= scalingFactors[i];
	}
	for (int i = 0; i < numberOfSplines; i++)
	{
		splineSumsNew[i] = max(0.0, splineSums[i] - sumOldPerBin[i] + sumNewPerBin[i]);
	}

	for (int i = 0; i < numberOfStandardParameters; i++)
	{
		sum += uR[i] * splineSumsNew[i];
	}
	for (int i = 0; i < numberOfSpecialParameters; i++)
	{
		sum += uR[N_PARAM - (numberOfSpecialParameters - i)] * (bcFactors[i][0] * splineSumsNew[numberOfSplines - 3] + bcFactors[i][1] * splineSumsNew[numberOfSplines - 2] + bcFactors[i][2] * splineSumsNew[numberOfSplines - 1] + bcFactors[i][3]);
	}

	exponentNew = sum;
	wfNew = exp(exponentNew + phiR);
}

double PBosonsBulk::CalculateWFQuotient(vector<vector<double> >& R, vector<double>& uR, vector<double>& uI, double phiR, double phiI, int changedParticleIndex, vector<double>& oldPosition)
{
	CalculateWFChange(R, uR, uI, phiR, phiI, changedParticleIndex, oldPosition);
	double wfQuotient = exp(2.0 * (exponentNew - exponent));

	//cout << "qu=" << wfQuotient << endl;

	return wfQuotient;
}

void PBosonsBulk::AcceptMove()
{
	wf = wfNew;
	exponent = exponentNew;
	splineSums = splineSumsNew;
}

void PBosonsBulk::FillCorrelatedSamplingData(ICorrelatedSamplingData* data)
{
	CSDataBulkSplines* d = dynamic_cast<CSDataBulkSplines*>(data);
	d->splineSumsD = this->splineSumsD;
	d->splineSumsD2 = this->splineSumsD2;
	d->otherLocalOperators = this->otherLocalOperators;
	d->grBins = this->grBins;
}
