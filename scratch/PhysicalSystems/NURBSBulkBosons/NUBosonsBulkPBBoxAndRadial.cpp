#include "NUBosonsBulkPBBoxAndRadial.h"

#include "../SplineFactory.h"

namespace PhysicalSystems
{

NUBosonsBulkPBBoxAndRadial::NUBosonsBulkPBBoxAndRadial(vector<double>& params, string configDirectory) :
		IPhysicalSystem(params, configDirectory)
{
	this->USE_NORMALIZATION_AND_PHASE = true;
	this->USE_NIC = true;
	this->USE_MOVE_COM_TO_ZERO = false;

	numberOfSplines = 0;
	numberOfSplinesRad = 0;

	halfLength = 0;
	maxDistanceRad = 0;
	numberOfSpecialParameters = 0;
	numberOfStandardParameters = 0;
	numberOfStandardParametersRad = 0;

	numOfOtherLocalOperators = 0;
	grBinCount = 0;
	grBinStartIndex = 0;
	grMaxDistance = 0;
	grNodePointSpacing = 0;
	numOfkValues = 0;

	exponentNew = 0;
	changedParticleIndex = 0;
}

void NUBosonsBulkPBBoxAndRadial::SetNodes(vector<double> n)
{
	this->nodes = n;
	int periodicNodeCount = 3;
	for (int i = 0; i < periodicNodeCount; i++)
	{
		this->nodes.insert(this->nodes.begin(), -n[i + 1]);
		this->nodes.push_back(2.0 * n[n.size() - 1] - n[n.size() - 2 - i]);
	}

	//this->nodesRad = this->nodes;
	//for (int i = 0; i < 0; i++)
	//{
	//	this->nodesRad.pop_back();
	//}
	//n.clear();
	//for (int i = 0; i < 31; i++)
	//{
	//	n.push_back(i * 0.05);
	//}
	this->nodesRad = n;
	for (int i = 0; i < periodicNodeCount; i++)
	{
		this->nodesRad.insert(this->nodesRad.begin(), -n[i + 1]);
		this->nodesRad.push_back(2.0 * n[n.size() - 1] - n[n.size() - 2 - i]);
	}
}

void NUBosonsBulkPBBoxAndRadial::SetGrBinCount(double n)
{
	this->grBinCount = n;
}

void NUBosonsBulkPBBoxAndRadial::InitSystem()
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

	//INFO: BC
	numberOfStandardParametersRad = N_PARAM / 2;
	numberOfStandardParameters = N_PARAM / 2;
	//numberOfStandardParametersRad = 30;
	//numberOfStandardParameters = N_PARAM - numberOfStandardParametersRad;
	numberOfSplinesRad = numberOfStandardParametersRad + 3;
	numberOfSplines = numberOfStandardParameters + 3;
	halfLength = LBOX / 2.0;
	//maxDistance = nodes[nodes.size() - 4];
	//maxDistance = halfLength;
	maxDistanceRad = nodesRad[nodesRad.size() - 4];
	//maxDistanceRad = 1.5;//nodesRad[nodesRad.size() - 4];
	//halfLength = maxDistanceRad;
	//halfLength = nodes[nodes.size() - 4];

	if (nodes.empty())
	{
		for (int i = -3; i < N_PARAM + 3; i++)
		{
			nodes.push_back((i * halfLength) / ((double) N_PARAM - 1));
		}
	}
	splineWeightsRad = SplineFactory::GetWeights(nodesRad);
	splineWeights = SplineFactory::GetWeights(nodes);

	//cut off
	//SplineFactory::SetBoundaryConditions1_1(nodes, bcFactors);
	//SplineFactory::SetBoundaryConditions1_2(nodes, bcFactors);
	//SplineFactory::SetBoundaryConditions2_1(nodes, bcFactors);
	//SplineFactory::SetBoundaryConditions2_2(nodes, bcFactors);
	//numberOfSpecialParameters = bcFactors.size();
	//numberOfStandardParameters = N_PARAM - numberOfSpecialParameters;
	//numberOfStandardParametersRad = (N_PARAM - 1) / 2;
	//numberOfStandardParametersRad = (N_PARAM -10)/ 2;
	//numberOfStandardParameters = N_PARAM - numberOfStandardParametersRad;

	wf = 0.0;
	exponent = 0.0;
	InitVector(splineSums, numberOfSplines, 0.0);
	InitVector(splineSumsD, numberOfSplines, N, DIM, 0.0);
	InitVector(splineSumsD2, numberOfSplines, N, 0.0);
	InitVector(splineSumsRad, numberOfSplinesRad, 0.0);
	InitVector(splineSumsDRad, numberOfSplinesRad, N, DIM, 0.0);
	InitVector(splineSumsD2Rad, numberOfSplinesRad, N, 0.0);

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
	InitVector(splineSumsNewRad, numberOfSplinesRad, 0.0);
	InitVector(sumOldPerBinRad, numberOfSplinesRad, 0.0);
	InitVector(sumNewPerBinRad, numberOfSplinesRad, 0.0);

	InitVector(grBins, grBinCount, 0.0);
	grMaxDistance = halfLength;
	grNodePointSpacing = grMaxDistance / (double) grBinCount;

	InitVector(grBinVolumes, grBinCount, 0.0);
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
	//pairDistribution.InitGrid(0.0, halfLength, 0.05);
	pairDistribution.InitGrid(0.0, halfLength, 0.02);
	pairDistribution.InitScaling();
	pairDistribution.InitObservables( { "g_2(r_ij)" });

	additionalObservables.Add(&pairDistribution);
}

void NUBosonsBulkPBBoxAndRadial::RefreshLocalOperators()
{
	//INFO: BC
	for (int i = 0; i < numberOfStandardParametersRad; i++)
	{
		this->localOperators[i] = splineSumsRad[i + 1];
	}
	this->localOperators[1] += splineSumsRad[0];
	this->localOperators[numberOfStandardParametersRad - 1] += splineSumsRad[numberOfSplinesRad - 2] / (-2.0);
	this->localOperators[numberOfStandardParametersRad - 1] += splineSumsRad[numberOfSplinesRad - 1];
	for (int i = 0; i < numberOfStandardParameters; i++)
	{
		this->localOperators[i + numberOfStandardParametersRad] = splineSums[i + 1];
	}
	this->localOperators[1 + numberOfStandardParametersRad] += splineSums[0];
	this->localOperators[N_PARAM - 1] += splineSums[numberOfSplines - 2];
	this->localOperators[N_PARAM - 1] += splineSums[numberOfSplines - 1];
}

void NUBosonsBulkPBBoxAndRadial::CalculateLocalOperators(vector<vector<double> >& R)
{
	int bin;
	double rni;
	double rni2;
	double rni3;
	double rnia;
	double rnia2;
	double rnia3;
	vector<double> vecrni(DIM);

	ClearVector(splineSums);
	ClearVector(splineSumsRad);

	for (int n = 0; n < N; n++)
	{
		for (int i = 0; i < n; i++)
		{
			rni = VectorDisplacementNIC_DIM(R[n], R[i], vecrni);
			if (rni < maxDistanceRad)
			{
				auto interval = lower_bound(nodesRad.begin(), nodesRad.end(), rni);
				bin = interval - nodesRad.begin() - 1;
				rni2 = rni * rni;
				rni3 = rni2 * rni;

				for (int p = 0; p < 4; p++)
				{
					splineSumsRad[bin - p] += splineWeightsRad[bin - p][p][0] + splineWeightsRad[bin - p][p][1] * rni + splineWeightsRad[bin - p][p][2] * rni2 + splineWeightsRad[bin - p][p][3] * rni3;
				}
			}
			//else
			{
				for (int a = 0; a < DIM; a++)
				{
					rnia = abs(vecrni[a]);
					auto interval = lower_bound(nodes.begin(), nodes.end(), rnia);
					bin = interval - nodes.begin() - 1;
					rnia2 = rnia * rnia;
					rnia3 = rnia2 * rnia;

					for (int p = 0; p < 4; p++)
					{
						splineSums[bin - p] += splineWeights[bin - p][p][0] + splineWeights[bin - p][p][1] * rnia + splineWeights[bin - p][p][2] * rnia2 + splineWeights[bin - p][p][3] * rnia3;
					}
				}
			}
		}
	}

	RefreshLocalOperators();
}

void NUBosonsBulkPBBoxAndRadial::CalculateOtherLocalOperators(vector<vector<double> >& R)
{
	int s = 4;
	vector<double> vecrni(DIM);
	vector<double> evecrni(DIM);
	vector<double> tmp1(s);
	vector<double> tmp2(s);
	int bin;
	double rni;
	double rni2;
	double rnia;
	double rnia2;

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
	ClearVector(splineSumsDRad);
	ClearVector(splineSumsD2Rad);
	ClearVector(grBins);
	ClearVector(otherLocalOperators);

	for (int n = 0; n < N; n++)
	{
		//external potential energy
		potentialExtern += GetExternalPotential(R[n]);

		for (int i = 0; i < N; i++)
		{
			rni = VectorDisplacementNIC_DIM(R[n], R[i], vecrni);

			//binning for g(r)
			if (i < n && rni < grMaxDistance)
			{
				grBinInterval = rni / grNodePointSpacing;
				grBin = int(grBinInterval);
				//grBins[grBin] += 1.0 / (double)pow(grBin + 1, 2);
				grBins[grBin] += 1.0 / grBinVolumes[grBin];
			}

			//internal potential energy
			if (rni < maxDistanceRad)
			{
				if (i < n)
				{
					//gauss potential
					double rnia = rni / a;
					potentialIntern += b * exp(-(rnia * rnia) / 2.0);

					//step potential
					//if (rni < a)
					//{
					//	potentialIntern += b;
					//}

					//double rniAbsorption = rni - 4.5;
					//double absorptionFactor = 0.0;//-3e-6;
					//if (rniAbsorption > 0)
					//{
					//	potentialInternComplex += absorptionFactor * rniAbsorption * rniAbsorption;
					//}
				}
			}

			//kinetic energy
			if (rni < maxDistanceRad)
			{
				if (i != n) //TODO: also include in (i < n) branch -> then only loop over i < n in "for (int i = 0; i < N; i++)"
				{
					auto interval = lower_bound(nodesRad.begin(), nodesRad.end(), rni);
					bin = interval - nodesRad.begin() - 1;
					rni2 = rni * rni;

					for (int p = 0; p < s; p++)
					{
						//first derivative
						tmp1[3 - p] = splineWeightsRad[bin - p][p][1] + 2.0 * splineWeightsRad[bin - p][p][2] * rni + 3.0 * splineWeightsRad[bin - p][p][3] * rni2;
						//second derivative
						tmp2[3 - p] = 2.0 * splineWeightsRad[bin - p][p][2] + 6.0 * splineWeightsRad[bin - p][p][3] * rni;
					}

					for (int a = 0; a < DIM; a++)
					{
						evecrni[a] = vecrni[a] / rni;
					}
					for (int a = 0; a < DIM; a++)
					{
						for (int b = 0; b < s; b++)
						{
							splineSumsDRad[bin - b][n][a] += tmp1[3 - b] * evecrni[a];
						}
					}
					double secondDerivativeFactor = DIM - 1.0; //INFO: at least correct for 2D and 3D (and 1D?)
					for (int b = 0; b < s; b++)
					{
						splineSumsD2Rad[bin - b][n] += tmp2[3 - b] + secondDerivativeFactor / rni * tmp1[3 - b];
					}
				}
			}
			//else
			{
				for (int a = 0; a < DIM; a++)
				{
					rnia = abs(vecrni[a]);
					//kinetic energy
					//TODO: maybe it is faster to calculate all localOperators[i][n] before and then loop over this precalculated values
					//		otherwise values are calculated multiple times
					if (i != n) //TODO: also include in (i < n) branch -> then only loop over i < n in "for (int i = 0; i < N; i++)"
					{
						auto interval = lower_bound(nodes.begin(), nodes.end(), rnia);
						bin = interval - nodes.begin() - 1;
						rnia2 = rnia * rnia;

						for (int p = 0; p < s; p++)
						{
							//first derivative
							tmp1[3 - p] = splineWeights[bin - p][p][1] + 2.0 * splineWeights[bin - p][p][2] * rnia + 3.0 * splineWeights[bin - p][p][3] * rnia2;
							//second derivative
							tmp2[3 - p] = 2.0 * splineWeights[bin - p][p][2] + 6.0 * splineWeights[bin - p][p][3] * rnia;
						}

						//for (int a = 0; a < DIM; a++)
						{
							for (int b = 0; b < s; b++)
							{
								int sign = vecrni[a] < 0 ? -1.0 : 1.0;
								splineSumsD[bin - b][n][a] += tmp1[3 - b] * sign;
							}
						}
						for (int b = 0; b < s; b++)
						{
							splineSumsD2[bin - b][n] += tmp2[3 - b];
						}
					}
				}
			}
		}
	}
	otherLocalOperators[0] = potentialIntern;
	otherLocalOperators[1] = potentialInternComplex;
	otherLocalOperators[2] = potentialExtern;
	otherLocalOperators[3] = potentialExternComplex;
}

vector<double> NUBosonsBulkPBBoxAndRadial::GetCenterOfMass(vector<vector<double> >& R)
{
	vector<double> com = { 0.0, 0.0, 0.0 };
	return com;
}

double NUBosonsBulkPBBoxAndRadial::GetExternalPotential(vector<double>& r)
{
	return 0;
}

void NUBosonsBulkPBBoxAndRadial::CalculateExpectationValues(vector<double>& O, vector<vector<vector<double> > >& sD, vector<vector<double> >& sD2, vector<vector<vector<double> > >& sDRad, vector<vector<double> >& sD2Rad, vector<double>& otherO, vector<double>& gr, vector<double>& uR, vector<double>& uI, double phiR, double phiI)
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

		//INFO: BC
		for (int k = 0; k < numberOfStandardParametersRad; k++)
		{
			for (int a = 0; a < DIM; a++)
			{
				vecKineticSumR1[a] += uR[k] * sDRad[k + 1][n][a];
				vecKineticSumI1[a] += uI[k] * sDRad[k + 1][n][a];
			}
			kineticSumR2 += uR[k] * sD2Rad[k + 1][n];
			kineticSumI2 += uI[k] * sD2Rad[k + 1][n];
		}
		for (int a = 0; a < DIM; a++)
		{
			vecKineticSumR1[a] += uR[1] * sDRad[0][n][a];
			vecKineticSumI1[a] += uI[1] * sDRad[0][n][a];
		}
		kineticSumR2 += uR[1] * sD2Rad[0][n];
		kineticSumI2 += uI[1] * sD2Rad[0][n];
		for (int a = 0; a < DIM; a++)
		{
			vecKineticSumR1[a] += uR[numberOfStandardParametersRad - 1] * sDRad[numberOfSplinesRad - 2][n][a] / (-2.0);
			vecKineticSumI1[a] += uI[numberOfStandardParametersRad - 1] * sDRad[numberOfSplinesRad - 2][n][a] / (-2.0);
		}
		kineticSumR2 += uR[numberOfStandardParametersRad - 1] * sD2Rad[numberOfSplinesRad - 2][n] / (-2.0);
		kineticSumI2 += uI[numberOfStandardParametersRad - 1] * sD2Rad[numberOfSplinesRad - 2][n] / (-2.0);
		for (int a = 0; a < DIM; a++)
		{
			vecKineticSumR1[a] += uR[numberOfStandardParametersRad - 1] * sD[numberOfSplinesRad - 1][n][a];
			vecKineticSumI1[a] += uI[numberOfStandardParametersRad - 1] * sD[numberOfSplinesRad - 1][n][a];
		}
		kineticSumR2 += uR[numberOfStandardParametersRad - 1] * sD2Rad[numberOfSplinesRad - 1][n];
		kineticSumI2 += uI[numberOfStandardParametersRad - 1] * sD2Rad[numberOfSplinesRad - 1][n];
		for (int k = 0; k < numberOfStandardParameters; k++)
		{
			for (int a = 0; a < DIM; a++)
			{
				vecKineticSumR1[a] += uR[k + numberOfStandardParametersRad] * sD[k + 1][n][a];
				vecKineticSumI1[a] += uI[k + numberOfStandardParametersRad] * sD[k + 1][n][a];
			}
			kineticSumR2 += uR[k + numberOfStandardParametersRad] * sD2[k + 1][n];
			kineticSumI2 += uI[k + numberOfStandardParametersRad] * sD2[k + 1][n];
		}
		for (int a = 0; a < DIM; a++)
		{
			vecKineticSumR1[a] += uR[1 + numberOfStandardParametersRad] * sD[0][n][a];
			vecKineticSumI1[a] += uI[1 + numberOfStandardParametersRad] * sD[0][n][a];
		}
		kineticSumR2 += uR[1 + numberOfStandardParametersRad] * sD2[0][n];
		kineticSumI2 += uI[1 + numberOfStandardParametersRad] * sD2[0][n];
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

void NUBosonsBulkPBBoxAndRadial::CalculateExpectationValues(vector<vector<double> >& R, vector<double>& uR, vector<double>& uI, double phiR, double phiI)
{
	RefreshLocalOperators();
	CalculateOtherLocalOperators(R);
	CalculateExpectationValues(this->localOperators, this->splineSumsD, this->splineSumsD2, this->splineSumsDRad, this->splineSumsD2Rad, this->otherLocalOperators, this->grBins, uR, uI, phiR, phiI);
}

void NUBosonsBulkPBBoxAndRadial::CalculateExpectationValues(ICorrelatedSamplingData* sample, vector<double>& uR, vector<double>& uI, double phiR, double phiI)
{
	CSDataBulkSplinesBR* s = dynamic_cast<CSDataBulkSplinesBR*>(sample);
	CalculateExpectationValues(s->localOperators, s->splineSumsD, s->splineSumsD2, s->splineSumsDRad, s->splineSumsD2Rad, s->otherLocalOperators, s->grBins, uR, uI, phiR, phiI);
}

void NUBosonsBulkPBBoxAndRadial::CalculateAdditionalSystemProperties(vector<vector<double> >& R, vector<double>& uR, vector<double>& uI, double phiR, double phiI)
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
}

void NUBosonsBulkPBBoxAndRadial::CalculateWavefunction(vector<double>& O, vector<double>& uR, vector<double>& uI, double phiR, double phiI)
{
	double sum = 0;

	for (int i = 0; i < N_PARAM; i++)
	{
		sum += uR[i] * O[i];
	}

	exponent = sum;
	wf = exp(exponent + phiR);
}

void NUBosonsBulkPBBoxAndRadial::CalculateWavefunction(vector<vector<double> >& R, vector<double>& uR, vector<double>& uI, double phiR, double phiI)
{
	CalculateLocalOperators(R);
	CalculateWavefunction(this->localOperators, uR, uI, phiR, phiI);
}

void NUBosonsBulkPBBoxAndRadial::CalculateWavefunction(ICorrelatedSamplingData* sample, vector<double>& uR, vector<double>& uI, double phiR, double phiI)
{
	CSDataBulkSplinesBR* s = dynamic_cast<CSDataBulkSplinesBR*>(sample);
	CalculateWavefunction(s->localOperators, uR, uI, phiR, phiI);
}

void NUBosonsBulkPBBoxAndRadial::CalculateWFChange(vector<vector<double> >& R, vector<double>& uR, vector<double>& uI, double phiR, double phiI, int changedParticleIndex, vector<double>& oldPosition)
{
	double sum = 0;
	int bin;
	double rni;
	double rni2;
	double rni3;
	double rnia;
	double rnia2;
	double rnia3;
	vector<double> vecrni(DIM);

	ClearVector(sumOldPerBin);
	ClearVector(sumNewPerBin);
	ClearVector(sumOldPerBinRad);
	ClearVector(sumNewPerBinRad);
	for (int i = 0; i < N; i++)
	{
		if (i != changedParticleIndex)
		{
			rni = VectorDisplacementNIC_DIM(R[i], oldPosition, vecrni);
			if (rni < maxDistanceRad)
			{
				auto interval = lower_bound(nodesRad.begin(), nodesRad.end(), rni);
				bin = interval - nodesRad.begin() - 1;
				rni2 = rni * rni;
				rni3 = rni2 * rni;

				for (int p = 0; p < 4; p++)
				{
					sumOldPerBinRad[bin - p] += splineWeightsRad[bin - p][p][0] + splineWeightsRad[bin - p][p][1] * rni + splineWeightsRad[bin - p][p][2] * rni2 + splineWeightsRad[bin - p][p][3] * rni3;
				}
			}
			//else
			{
				for (int a = 0; a < DIM; a++)
				{
					rnia = abs(vecrni[a]);
					auto interval = lower_bound(nodes.begin(), nodes.end(), rnia);
					bin = interval - nodes.begin() - 1;
					rnia2 = rnia * rnia;
					rnia3 = rnia2 * rnia;

					for (int p = 0; p < 4; p++)
					{
						sumOldPerBin[bin - p] += splineWeights[bin - p][p][0] + splineWeights[bin - p][p][1] * rnia + splineWeights[bin - p][p][2] * rnia2 + splineWeights[bin - p][p][3] * rnia3;
					}
				}
			}
			rni = VectorDisplacementNIC_DIM(R[i], R[changedParticleIndex], vecrni);
			if (rni < maxDistanceRad)
			{
				auto interval = lower_bound(nodesRad.begin(), nodesRad.end(), rni);
				bin = interval - nodesRad.begin() - 1;
				rni2 = rni * rni;
				rni3 = rni2 * rni;

				for (int p = 0; p < 4; p++)
				{
					sumNewPerBinRad[bin - p] += splineWeightsRad[bin - p][p][0] + splineWeightsRad[bin - p][p][1] * rni + splineWeightsRad[bin - p][p][2] * rni2 + splineWeightsRad[bin - p][p][3] * rni3;
				}
			}
			//else
			{
				for (int a = 0; a < DIM; a++)
				{
					rnia = abs(vecrni[a]);
					auto interval = lower_bound(nodes.begin(), nodes.end(), rnia);
					bin = interval - nodes.begin() - 1;
					rnia2 = rnia * rnia;
					rnia3 = rnia2 * rnia;

					for (int p = 0; p < 4; p++)
					{
						sumNewPerBin[bin - p] += splineWeights[bin - p][p][0] + splineWeights[bin - p][p][1] * rnia + splineWeights[bin - p][p][2] * rnia2 + splineWeights[bin - p][p][3] * rnia3;
					}
				}
			}
		}
	}

	for (int i = 0; i < numberOfSplinesRad; i++)
	{
		splineSumsNewRad[i] = max(0.0, splineSumsRad[i] - sumOldPerBinRad[i] + sumNewPerBinRad[i]);
	}
	for (int i = 0; i < numberOfSplines; i++)
	{
		splineSumsNew[i] = max(0.0, splineSums[i] - sumOldPerBin[i] + sumNewPerBin[i]);
	}

	//INFO: BC
	for (int i = 0; i < numberOfStandardParametersRad; i++)
	{
		sum += uR[i] * splineSumsNewRad[i + 1];
	}
	sum += uR[1] * splineSumsNewRad[0];
	sum += uR[numberOfStandardParametersRad - 1] * splineSumsNewRad[numberOfSplinesRad - 2] / (-2.0);
	sum += uR[numberOfStandardParametersRad - 1] * splineSumsNewRad[numberOfSplinesRad - 1];
	for (int i = 0; i < numberOfStandardParameters; i++)
	{
		sum += uR[i + numberOfStandardParametersRad] * splineSumsNew[i + 1];
	}
	sum += uR[1 + numberOfStandardParametersRad] * splineSumsNew[0];
	sum += uR[N_PARAM - 1] * splineSumsNew[numberOfSplines - 2];
	sum += uR[N_PARAM - 1] * splineSumsNew[numberOfSplines - 1];

	exponentNew = sum;
	wfNew = exp(exponentNew + phiR);
}

double NUBosonsBulkPBBoxAndRadial::CalculateWFQuotient(vector<vector<double> >& R, vector<double>& uR, vector<double>& uI, double phiR, double phiI, int changedParticleIndex, vector<double>& oldPosition)
{
	CalculateWFChange(R, uR, uI, phiR, phiI, changedParticleIndex, oldPosition);
	double wfQuotient = exp(2.0 * (exponentNew - exponent));

	//cout << "qu=" << wfQuotient << " (" << exponentNew << " - " << exponent << ")" << endl;

	return wfQuotient;
}

void NUBosonsBulkPBBoxAndRadial::AcceptMove()
{
	wf = wfNew;
	exponent = exponentNew;
	splineSums = splineSumsNew;
	splineSumsRad = splineSumsNewRad;
}

void NUBosonsBulkPBBoxAndRadial::InitCorrelatedSamplingData(vector<ICorrelatedSamplingData*>& data)
{
	for (unsigned int i = 0; i < data.size(); i++)
	{
		data[i] = new CSDataBulkSplinesBR();
	}
}

void NUBosonsBulkPBBoxAndRadial::FillCorrelatedSamplingData(ICorrelatedSamplingData* data)
{
	CSDataBulkSplinesBR* d = dynamic_cast<CSDataBulkSplinesBR*>(data);
	d->splineSumsD = this->splineSumsD;
	d->splineSumsD2 = this->splineSumsD2;
	d->splineSumsDRad = this->splineSumsDRad;
	d->splineSumsD2Rad = this->splineSumsD2Rad;
	d->otherLocalOperators = this->otherLocalOperators;
	d->grBins = this->grBins;
}

}
