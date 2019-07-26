#include "NUBosonsBulkPB.h"

#include "../SplineFactory.h"

namespace PhysicalSystems
{

NUBosonsBulkPB::NUBosonsBulkPB(vector<double>& params, string configDirectory) :
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

	outerSum = 0.0;
	outerSumNew = 0.0;
	exponentNew = 0;
	changedParticleIndex = 0;
}

void NUBosonsBulkPB::SetNodes(vector<double> n)
{
	this->nodes = n;
	int periodicNodeCount = 3;
	for (int i = 0; i < periodicNodeCount; i++)
	{
		this->nodes.insert(this->nodes.begin(), -n[i + 1]);
		this->nodes.push_back(2.0 * n[n.size() - 1] - n[n.size() - 2 - i]);
	}
}

void NUBosonsBulkPB::SetGrBinCount(double n)
{
	this->grBinCount = n;
}

void NUBosonsBulkPB::InitSystem()
{
	numOfOtherLocalOperators = 4; //INFO: potentialIntern, potentialInternComplex, potentialExtern, potentialExternComplex
	otherLocalOperators.resize(numOfOtherLocalOperators);
	if (grBinCount == 0)
	{
		grBinCount = 400;
	}
	grBinStartIndex = 3;
	numOfkValues = 300;
	numOfOtherExpectationValues = 3 + grBinCount;
	numOfAdditionalSystemProperties = numOfOtherExpectationValues + numOfkValues;
	//numOfAdditionalSystemProperties = numOfOtherExpectationValues + numOfkValues + numberOfSplines + numberOfSplines * N * DIM + numberOfSplines * N; //expectationValues + s(k) + splines + splinesD + splinesD2

	//INFO: (N_PARAM + 2) for imaginary time with BC1.1;
	//INFO: (N_PARAM + 1) for BC2.2
	//INFO: N_PARAM + 2 - 2 for no BC and real time;
	//INFO: +1 for ignoring first one
	//numberOfSplines = N_PARAM + 2;//BC 1.1
	numberOfSplines = N_PARAM + 3;//BC 1.2
	//numberOfSplines = N_PARAM + 2;//BC 2.1
	halfLength = LBOX / 2.0;
	//maxDistance = nodes[nodes.size() - 4];
	maxDistance = halfLength;

	//nodes.resize(numberOfSplines + 4);
	//double baseSpacing = pow(maxDistance, 1.0/1.5) / (double) (N_PARAM - 2.0);
	//nodes[0] = -baseSpacing * 2.0;
	//nodes[1] = -baseSpacing;
	//nodes[2] = -baseSpacing / 2.0;
	//nodes[3] = baseSpacing / 2.0;
	//nodes[4] = baseSpacing;
	//nodes[5] = baseSpacing * 2.0;
	//for (int i = 0; i < numberOfSplines + 4 - 3; i++)
	//{
	//	nodes[i + 6] = pow(i * baseSpacing, 1.5);
	//}
	if (nodes.empty())
	{
		for (int i = -3; i < N_PARAM + 3;  i++)
		{
			nodes.push_back((i * halfLength) / ((double)N_PARAM - 1));
		}
	}
	splineWeights = SplineFactory::GetWeights(nodes);

	//cut off
	//SplineFactory::SetBoundaryConditions1_1(nodes, bcFactors);
	//SplineFactory::SetBoundaryConditions1_2(nodes, bcFactors);
	//SplineFactory::SetBoundaryConditions2_1(nodes, bcFactors);
	//SplineFactory::SetBoundaryConditions2_2(nodes, bcFactors);
	numberOfSpecialParameters = bcFactors.size();
	numberOfStandardParameters = N_PARAM - numberOfSpecialParameters;

	wf = 0.0;
	exponent = 0.0;
	outerSum = 0.0;
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
		scalingFactors[i] = 1.0;
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
	grMaxDistance = halfLength;
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

	pairDistribution.name = "pairDistribution";
	pairDistribution.InitGrid(0.0, this->nodes[this->nodes.size() - 4], 0.05);
	pairDistribution.InitScaling();
	pairDistribution.InitObservables({"g_2(r_ij)"});

	additionalObservables.Add(&pairDistribution);
}

void NUBosonsBulkPB::RefreshLocalOperators()
{
	for (int i = 0; i < numberOfStandardParameters; i++)
	{
		this->localOperators[i] = splineSums[i + 1];
	}
	this->localOperators[1] += splineSums[0];
	this->localOperators[numberOfStandardParameters - 1] += splineSums[numberOfSplines - 2];
	this->localOperators[numberOfStandardParameters - 1] += splineSums[numberOfSplines - 1];
}

void NUBosonsBulkPB::CalculateLocalOperators(vector<vector<double> >& R)
{
	int bin;
	double rni;
	double rni2;
	double rni3;
	vector<double> vecrni(DIM);

	outerSum = 0.0;
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
				rni2 = rni * rni;
				rni3 = rni2 * rni;

				for (int p = 0; p < 4; p++)
				{
					splineSums[bin - p] += splineWeights[bin - p][p][0] + splineWeights[bin - p][p][1] * rni + splineWeights[bin - p][p][2] * rni2 + splineWeights[bin - p][p][3] * rni3;
				}
			}
			else
			{
				outerSum += 1.0;
			}
		}
	}
	for (int i = 0; i < numberOfSplines; i++)
	{
		splineSums[i] *= scalingFactors[i];
	}
	RefreshLocalOperators();
}

void NUBosonsBulkPB::CalculateOtherLocalOperators(vector<vector<double> >& R)
{
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
			rni = VectorDisplacementNIC_DIM(R[n], R[i], vecrni);
			if (rni < maxDistance)
			{
				//internal potential energy
				if (i < n)
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
				if (i != n) //TODO: also include in (i < n) branch -> then only loop over i < n in "for (int i = 0; i < N; i++)"
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
						splineSumsD2[bin - b][n] += tmp2[3 - b] + secondDerivativeFactor/ rni * tmp1[3 - b];
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

vector<double> NUBosonsBulkPB::GetCenterOfMass(vector<vector<double> >& R)
{
	vector<double> com = { 0.0, 0.0, 0.0 };
	return com;
}

double NUBosonsBulkPB::GetExternalPotential(vector<double>& r)
{
	return 0;
}

void NUBosonsBulkPB::CalculateExpectationValues(vector<double>& O, vector<vector<vector<double> > >& sD, vector<vector<double> >& sD2, vector<double>& otherO, vector<double>& gr, vector<double>& uR, vector<double>& uI, double phiR, double phiI)
{
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
				vecKineticSumR1[a] += uR[k] * sD[k + 1][n][a];
				vecKineticSumI1[a] += uI[k] * sD[k + 1][n][a];
			}
			kineticSumR2 += uR[k] * sD2[k + 1][n];
			kineticSumI2 += uI[k] * sD2[k + 1][n];
		}
		for (int a = 0; a < DIM; a++)
		{
			vecKineticSumR1[a] += uR[1] * sD[0][n][a];
			vecKineticSumI1[a] += uI[1] * sD[0][n][a];
		}
		kineticSumR2 += uR[1] * sD2[0][n];
		kineticSumI2 += uI[1] * sD2[0][n];
		for (int a = 0; a < DIM; a++)
		{
			vecKineticSumR1[a] += uR[N_PARAM - 1] * sD[numberOfSplines - 2][n][a];
			vecKineticSumI1[a] += uI[N_PARAM - 1] * sD[numberOfSplines - 2][n][a];
		}
		kineticSumR2 += uR[N_PARAM - 1] * sD2[numberOfSplines - 2][n];
		kineticSumI2 += uI[N_PARAM - 1] * sD2[numberOfSplines - 2][n];
		for (int a = 0; a < DIM; a++)
		{
			vecKineticSumR1[a] += uR[N_PARAM - 1] * sD[numberOfSplines - 1][n][a];
			vecKineticSumI1[a] += uI[N_PARAM - 1] * sD[numberOfSplines - 1][n][a];
		}
		kineticSumR2 += uR[N_PARAM - 1] * sD2[numberOfSplines - 1][n];
		kineticSumI2 += uI[N_PARAM - 1] * sD2[numberOfSplines - 1][n];

		kineticSumR1I1 += 2.0 * VectorDotProduct_DIM(vecKineticSumR1, vecKineticSumI1);
		kineticSumR1 += VectorNorm2_DIM(vecKineticSumR1);
		kineticSumI1 += VectorNorm2_DIM(vecKineticSumI1);
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

void NUBosonsBulkPB::CalculateExpectationValues(vector<vector<double> >& R, vector<double>& uR, vector<double>& uI, double phiR, double phiI)
{
	RefreshLocalOperators();
	CalculateOtherLocalOperators(R);
	CalculateExpectationValues(this->localOperators, this->splineSumsD, this->splineSumsD2, this->otherLocalOperators, this->grBins, uR, uI, phiR, phiI);
}

void NUBosonsBulkPB::CalculateExpectationValues(ICorrelatedSamplingData* sample, vector<double>& uR, vector<double>& uI, double phiR, double phiI)
{
	CSDataBulkSplines* s = dynamic_cast<CSDataBulkSplines*>(sample);
	CalculateExpectationValues(s->localOperators, s->splineSumsD, s->splineSumsD2, s->otherLocalOperators, s->grBins, uR, uI, phiR, phiI);
}

void NUBosonsBulkPB::CalculateAdditionalSystemProperties(vector<vector<double> >& R, vector<double>& uR, vector<double>& uI, double phiR, double phiI)
{
	additionalObservables.ClearValues();
	double rni;
	vector<double> vecrni(DIM);

	//vector<double> sumSCos;
	//vector<double> sumSSin;
	//vector<double> sk;
	//sumSCos.resize(numOfkValues);
	//ClearVector(sumSCos);
	//sumSSin.resize(numOfkValues);
	//ClearVector(sumSSin);
	//sk.resize(numOfkValues);
	//ClearVector(sk);
	//for (int i = 0; i < N; i++)
	//{
	//	for (int k = 0; k < numOfkValues; k++)
	//	{
	//		for (unsigned int kn = 0; kn < kValues[k].size(); kn++)
	//		{
	//			//sumSCos[k] += cos(kValues[k][kn][0] * R[i][0] + kValues[k][kn][1] * R[i][1] + kValues[k][kn][2] * R[i][2]);
	//			//sumSSin[k] += sin(kValues[k][kn][0] * R[i][0] + kValues[k][kn][1] * R[i][1] + kValues[k][kn][2] * R[i][2]);
	//		}
	//	}
	//}
	//for (int k = 0; k < numOfkValues; k++)
	//{
	//	sk[k] = (sumSCos[k] * sumSCos[k] + sumSSin[k] * sumSSin[k]) / ((double) (N * kValues[k].size()));
	//	additionalSystemProperties[numOfOtherExpectationValues + k] = sk[k];
	//	//index++;
	//}

	//pairDistribution
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < i; j++)
		{
			rni = VectorDisplacementNIC_DIM(R[i], R[j], vecrni);
			if (rni < pairDistribution.gridMax)
			{
				pairDistribution.AddToHistogram(0, rni, 1.0);
			}
		}
	}
}

void NUBosonsBulkPB::CalculateWavefunction(vector<double>& O, vector<double>& uR, vector<double>& uI, double phiR, double phiI)
{
	double sum = 0;

	for (int i = 0; i < N_PARAM; i++)
	{
		sum += uR[i] * O[i];
	}
	//TODO: does not work for ICorrelatedSamplingData
	sum += uR.back() * outerSum;

	exponent = sum;
	wf = exp(exponent + phiR);
}

void NUBosonsBulkPB::CalculateWavefunction(vector<vector<double> >& R, vector<double>& uR, vector<double>& uI, double phiR, double phiI)
{
	CalculateLocalOperators(R);
	CalculateWavefunction(this->localOperators, uR, uI, phiR, phiI);
}

void NUBosonsBulkPB::CalculateWavefunction(ICorrelatedSamplingData* sample, vector<double>& uR, vector<double>& uI, double phiR, double phiI)
{
	CSDataBulkSplines* s = dynamic_cast<CSDataBulkSplines*>(sample);
	CalculateWavefunction(s->localOperators, uR, uI, phiR, phiI);
}

void NUBosonsBulkPB::CalculateWFChange(vector<vector<double> >& R, vector<double>& uR, vector<double>& uI, double phiR, double phiI, int changedParticleIndex, vector<double>& oldPosition)
{
	double sum = 0;
	int bin;
	double rni;
	double rni2;
	double rni3;
	vector<double> vecrni(DIM);

	double oldOuterSum = 0.0;
	double newOuterSum = 0.0;
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
				rni2 = rni * rni;
				rni3 = rni2 * rni;

				for (int p = 0; p < 4; p++)
				{
					sumOldPerBin[bin - p] += splineWeights[bin - p][p][0] + splineWeights[bin - p][p][1] * rni + splineWeights[bin - p][p][2] * rni2 + splineWeights[bin - p][p][3] * rni3;
				}
			}
			else
			{
				oldOuterSum += 1.0;
			}
			rni = VectorDisplacementNIC_DIM(R[i], R[changedParticleIndex], vecrni);
			if (rni < maxDistance)
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
			else
			{
				newOuterSum += 1.0;
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
	outerSumNew = max(0.0, outerSum - oldOuterSum + newOuterSum);

	for (int i = 0; i < numberOfStandardParameters; i++)
	{
		sum += uR[i] * splineSumsNew[i + 1];
	}
	sum += uR[1] * splineSumsNew[0];
	sum += uR[N_PARAM - 1] * splineSumsNew[numberOfSplines - 2];
	sum += uR[N_PARAM - 1] * splineSumsNew[numberOfSplines - 1];
	sum += uR.back() * outerSumNew;

	exponentNew = sum;
	wfNew = exp(exponentNew + phiR);
}

double NUBosonsBulkPB::CalculateWFQuotient(vector<vector<double> >& R, vector<double>& uR, vector<double>& uI, double phiR, double phiI, int changedParticleIndex, vector<double>& oldPosition)
{
	CalculateWFChange(R, uR, uI, phiR, phiI, changedParticleIndex, oldPosition);
	double wfQuotient = exp(2.0 * (exponentNew - exponent));

	//cout << "qu=" << wfQuotient << " (" << exponentNew << " - " << exponent << ")" << endl;

	return wfQuotient;
}

void NUBosonsBulkPB::AcceptMove()
{
	wf = wfNew;
	exponent = exponentNew;
	splineSums = splineSumsNew;
	outerSum = outerSumNew;
}

void NUBosonsBulkPB::FillCorrelatedSamplingData(ICorrelatedSamplingData* data)
{
	CSDataBulkSplines* d = dynamic_cast<CSDataBulkSplines*>(data);
	d->splineSumsD = this->splineSumsD;
	d->splineSumsD2 = this->splineSumsD2;
	d->otherLocalOperators = this->otherLocalOperators;
	d->grBins = this->grBins;
}

}
