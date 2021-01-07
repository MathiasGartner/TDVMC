#include "Bosons1D.h"

#include "../SplineFactory.h"

namespace PhysicalSystems
{

Bosons1D::Bosons1D(vector<double>& params, string configDirectory) :
		IPhysicalSystem(params, configDirectory)
{
	this->USE_NORMALIZATION_AND_PHASE = true;
	this->USE_NIC = true;
	this->USE_MOVE_COM_TO_ZERO = false;

	usePreCalcSplineValues = false;

	numberOfSplines = 0;

	halfLength = 0;
	maxDistance = 0;
	numberOfSpecialParametersStart = 0;
	numberOfSpecialParametersEnd = 0;
	numberOfStandardParameters = 0;
	np1 = 0;
	np2 = 0;
	np3 = 0;

	numOfOtherLocalOperators = 0;
	numOfPairDistributionValues = 0;
	numOfkValues = 0;

	exponentNew = 0;
	changedParticleIndex = 0;

	particleStartIndexForReducedSampling = 0;
}

void Bosons1D::SetNodes(vector<double> n)
{
	this->nodes = n;
	int periodicNodeCount = 3;
	for (int i = 0; i < periodicNodeCount; i++)
	{
		this->nodes.insert(this->nodes.begin(), -n[i + 1]);
		this->nodes.push_back(2.0 * n[n.size() - 1] - n[n.size() - 2 - i]);
	}
}

void Bosons1D::SetPairDistributionBinCount(double n)
{
	this->numOfPairDistributionValues = n;
}

void Bosons1D::SetParticleStartIndexForReducedSampling(int n)
{
	this->particleStartIndexForReducedSampling = n;
}

double Bosons1D::CalculateOBDMKernel(vector<double>& r, vector<vector<double> >& R, vector<double>& uR, vector<double>& uI, double phiR, double phiI)
{
	double result = 0.0;

	vector<vector<double> > firstParticleR = {R[0], r};
	vector<double> sumParts;
	InitVector(sumParts, 2, 0.0);

	vector<double> tmpSplineSums;
	vector<double> tmpLocalOperators;

	for (int part = 0; part < 2; part++)
	{
		//according to Bosons1D::CalculateLocalOperators(vector<vector<double> >& R)
		int bin;
		double rni;
		double rni2;
		double rni3;
		vector<double> vecrni(DIM);

		InitVector(tmpSplineSums, this->splineSums.size(), 0.0);
		InitVector(tmpLocalOperators, this->localOperators.size(), 0.0);

		for (int n = particleStartIndexForReducedSampling; n < N; n++)
		{
			rni = VectorDisplacementNIC_DIM(R[n], firstParticleR[part], vecrni);
			if (rni <= maxDistance)
			{
				auto interval = lower_bound(nodes.begin(), nodes.end(), rni);
				bin = interval - nodes.begin() - 1;
				rni2 = rni * rni;
				rni3 = rni2 * rni;

				for (int p = 0; p < 4; p++)
				{
					tmpSplineSums[bin - p] += splineWeights[bin - p][p][0] + splineWeights[bin - p][p][1] * rni + splineWeights[bin - p][p][2] * rni2 + splineWeights[bin - p][p][3] * rni3;
				}
			}
			else
			{
				cout << "!!!!" << endl;
			}
		}

		//according to Bosons1D::RefreshLocalOperators()
		//INFO: BC
		for (int i = 0; i < np1; i++)
		{
			tmpLocalOperators[i] = (bcFactorsStart[i][0] * tmpSplineSums[0] + bcFactorsStart[i][1] * tmpSplineSums[1] + bcFactorsStart[i][2] * tmpSplineSums[2]);
		}
		int spIdx = 3;
		for (int i = np1; i < np2; i++)
		{
			tmpLocalOperators[i] = tmpSplineSums[spIdx];
			spIdx++;
		}
		int idx = 0;
		for (int i = np2; i < np3; i++)
		{
			tmpLocalOperators[i] = (bcFactorsEnd[idx][0] * tmpSplineSums[numberOfSplines - 3] + bcFactorsEnd[idx][1] * tmpSplineSums[numberOfSplines - 2] + bcFactorsEnd[idx][2] * tmpSplineSums[numberOfSplines - 1]);
			idx++;
		}

		//according to Bosons1D::CalculateWavefunction(vector<double>& O, vector<double>& uR, vector<double>& uI, double phiR, double phiI)
		for (int i = 0; i < N_PARAM; i++)
		{
			sumParts[part] += uR[i] * tmpLocalOperators[i];
		}
	}

	//result = exp(sumParts[0] + sumParts[1] + 2.0 * phiR);
	result = sumParts[0] + sumParts[1];
	//result = sumParts[0];

	return result;
}

void Bosons1D::InitSystem()
{
	numOfOtherLocalOperators = 4; //INFO: potentialIntern, potentialInternComplex, potentialExtern, potentialExternComplex
	otherLocalOperators.resize(numOfOtherLocalOperators);
	numOfkValues = 150;
	numOfOtherExpectationValues = 9;
	numOfAdditionalSystemProperties = numOfOtherExpectationValues;

	halfLength = LBOX / 2.0;

	if (nodes.empty())
	{
		//INFO: BC
		//this is for SetBoundaryConditions3_1D_OR_2 and SetBoundaryConditions3_1D_CO_2
		for (int i = -3; i < N_PARAM + 3; i++)
		{
			nodes.push_back((i * halfLength) / ((double) N_PARAM - 1));
		}
	}
	splineWeights = SplineFactory::GetWeights(nodes);

	//INFO: BC
	//numberOfSplines = N_PARAM + 1; //BC3.1D.CO.1
	//numberOfSplines = N_PARAM + 2; //BC3.1D.CO.2
	//numberOfSplines = N_PARAM + 3; //BC3.1D.CO.2 + BC3.1D.OR.1
	numberOfSplines = nodes.size() - 3 - 1; //3rd order splines -(d+1) (d...spline dimension)
	maxDistance = nodes[nodes.size() - 4];

	//cut off
	SplineFactory::SetBoundaryConditions3_1D_OR_2(nodes, bcFactorsStart, true);
	SplineFactory::SetBoundaryConditions3_1D_CO_2(nodes, bcFactorsEnd, maxDistance, true);
	numberOfSpecialParametersStart = bcFactorsStart.size();
	numberOfSpecialParametersEnd = bcFactorsEnd.size();
	numberOfStandardParameters = N_PARAM - numberOfSpecialParametersStart - numberOfSpecialParametersEnd;
	np1 = numberOfSpecialParametersStart;
	np2 = np1 + numberOfStandardParameters;
	np3 = np2 + numberOfSpecialParametersEnd;

	int requiredParams = numberOfSplines - (3 - numberOfSpecialParametersStart) - (3 - numberOfSpecialParametersEnd);
	if (N_PARAM != requiredParams)
	{
		cout << "!!! WRONG NUMBER OF PARAMETERS !!!" << endl;
		cout << "!!! need " << requiredParams << " parameters, got " << N_PARAM << " !!!" << endl;
		exit(0);
	}

	wf = 0.0;
	exponent = 0.0;
	InitVector(splineSums, numberOfSplines, 0.0);
	InitVector(splineSumsD, numberOfSplines, N, DIM, 0.0);
	InitVector(splineSumsD2, numberOfSplines, N, 0.0);

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
	pairDistribution.grid.name = "r_ij";
	//pairDistribution.InitGrid(0.0, halfLength, 0.05);
	pairDistribution.InitGrid(0.0, nodes[nodes.size() - 4], nodes[nodes.size() - 4] / numOfPairDistributionValues);
	pairDistribution.InitScaling();
	pairDistribution.InitObservables( { "g_2(r_ij)" });

	structureFactor.name = "structureFactor";
	structureFactor.grid.name = "k";
	structureFactor.InitGrid(kNorms);
	structureFactor.InitObservables( { "S(k)" });

	additionalObservables.Add(&pairDistribution);
	additionalObservables.Add(&structureFactor);

	//Precalculate spline values:
	if (usePreCalcSplineValues)
	{
		int preCalcBins = 2000;
		double rni, rni2, rni3;
		InitVector(nodeDiffs, nodes.size() - 1, 0.0);
		InitVector(binSizePerNode, nodes.size() - 1, 0.0);
		InitVector(binSizePerNode_R, nodes.size() - 1, 0.0);
		for (unsigned int n = 0; n < nodes.size() - 1; n++)
		{
			nodeDiffs[n] = nodes[n + 1] - nodes[n];
			binSizePerNode[n] = nodeDiffs[n] / preCalcBins;
			binSizePerNode_R[n] = 1.0 / binSizePerNode[n];
		}
		InitVector(preCalcSplineValues, splineWeights.size(), splineWeights[0].size(), preCalcBins, 0.0);
		for (unsigned int n = 3; n < nodes.size() - 4; n++)
		{
			rni = nodes[n] + binSizePerNode[n] / 2.0;
			for (int i = 0; i < preCalcBins; i++)
			{
				rni2 = rni * rni;
				rni3 = rni2 * rni;
				for (int p = 0; p < 4; p++)
				{
					preCalcSplineValues[n - p][p][i] = splineWeights[n - p][p][0] + splineWeights[n - p][p][1] * rni + splineWeights[n - p][p][2] * rni2 + splineWeights[n - p][p][3] * rni3;
				}
				rni += binSizePerNode[n];
			}
		}
	}
}

void Bosons1D::RefreshLocalOperators()
{
	//INFO: BC
	for (int i = 0; i < np1; i++)
	{
		this->localOperators[i] = (bcFactorsStart[i][0] * splineSums[0] + bcFactorsStart[i][1] * splineSums[1] + bcFactorsStart[i][2] * splineSums[2]);
	}
	int spIdx = 3;
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

void Bosons1D::CalculateLocalOperators(vector<vector<double> >& R)
{
	int bin;
	double rni;
	double rni2;
	double rni3;
	vector<double> vecrni(DIM);

	ClearVector(splineSums);

	for (int n = particleStartIndexForReducedSampling; n < N; n++)
	{
		for (int i = particleStartIndexForReducedSampling; i < n; i++)
		{
			rni = VectorDisplacementNIC_DIM(R[n], R[i], vecrni);
			if (rni <= maxDistance)
			{
				auto interval = lower_bound(nodes.begin(), nodes.end(), rni);
				bin = interval - nodes.begin() - 1;
				rni2 = rni * rni;
				rni3 = rni2 * rni;

				if (usePreCalcSplineValues)
				{
					int rniBin = (rni - nodes[bin]) * binSizePerNode_R[bin + 3];
					for (int p = 0; p < 4; p++)
					{
						splineSums[bin - p] += preCalcSplineValues[bin - p][p][rniBin];
					}
				}
				else
				{
					for (int p = 0; p < 4; p++)
					{
						splineSums[bin - p] += splineWeights[bin - p][p][0] + splineWeights[bin - p][p][1] * rni + splineWeights[bin - p][p][2] * rni2 + splineWeights[bin - p][p][3] * rni3;
					}
				}
			}
			else
			{
				cout << "!!!!" << endl;
			}
		}
	}
	//INFO: with preCalcsplineValues
	//double sum = 0;
	//double sum2 = 0;
	//for (unsigned int i = 0; i < splineSums.size(); i++)
	//{
	//	cout << (splineSums[i] - splineSumsPreCalc[i]) << "(" << splineSums[i] << " - " << splineSumsPreCalc[i] << ")" << endl;
	//	sum += (splineSums[i] - splineSumsPreCalc[i]);
	//	sum2 += sqrt((splineSums[i] - splineSumsPreCalc[i]) * (splineSums[i] - splineSumsPreCalc[i]));
	//}
	//cout << "sum=" << sum << ", sum2=" << sum2 << endl;
	RefreshLocalOperators();
}

void Bosons1D::CalculateOtherLocalOperators(vector<vector<double> >& R)
{
	int s = 4;
	vector<double> vecrni(DIM);
	vector<double> evecrni(DIM);
	vector<double> tmp1(s);
	vector<double> tmp2(s);
	int bin;
	double rni;
	double rni2;

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

	ClearVector(splineSumsD);
	ClearVector(splineSumsD2);
	ClearVector(otherLocalOperators);

	for (int n = 0; n < N; n++)
	{
		//external potential energy
		potentialExtern += GetExternalPotential(R[n]);

		for (int i = 0; i < N; i++)
		{
			rni = VectorDisplacementNIC_DIM(R[n], R[i], vecrni);

			if (rni <= maxDistance)
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
					bool outsideL2 = false;
					double firstDerivativeFactor = outsideL2 ? -1.0 : 1.0;
					for (int a = 0; a < DIM; a++)
					{
						for (int b = 0; b < s; b++)
						{
							splineSumsD[bin - b][n][a] += tmp1[3 - b] * evecrni[a] * firstDerivativeFactor;
						}
					}
					double secondDerivativeFactor = DIM - 1.0; //INFO: at least correct for 2D and 3D (and 1D?)
					for (int b = 0; b < s; b++)
					{
						splineSumsD2[bin - b][n] += tmp2[3 - b] + secondDerivativeFactor / rni * tmp1[3 - b];
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

vector<double> Bosons1D::GetCenterOfMass(vector<vector<double> >& R)
{
	vector<double> com = { 0.0, 0.0, 0.0 };
	return com;
}

double Bosons1D::GetExternalPotential(vector<double>& r)
{
	return 0;
}

void Bosons1D::CalculateExpectationValues(vector<double>& O, vector<vector<vector<double> > >& sD, vector<vector<double> >& sD2, vector<double>& otherO, vector<double>& gr, vector<double>& uR, vector<double>& uI, double phiR, double phiI)
{
	double tmp = 0;
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
		for (int i = 0; i < np1; i++)
		{
			for (int a = 0; a < DIM; a++)
			{
				tmp = (bcFactorsStart[i][0] * sD[0][n][a] + bcFactorsStart[i][1] * sD[1][n][a] + bcFactorsStart[i][2] * sD[2][n][a]);
				vecKineticSumR1[a] += uR[i] * tmp;
				vecKineticSumI1[a] += uI[i] * tmp;
			}
			tmp = (bcFactorsStart[i][0] * sD2[0][n] + bcFactorsStart[i][1] * sD2[1][n] + bcFactorsStart[i][2] * sD2[2][n]);
			kineticSumR2 += uR[i] * tmp;
			kineticSumI2 += uI[i] * tmp;
		}
		int spIdx = 3;
		for (int i = np1; i < np2; i++)
		{
			for (int a = 0; a < DIM; a++)
			{
				vecKineticSumR1[a] += uR[i] * sD[spIdx][n][a];
				vecKineticSumI1[a] += uI[i] * sD[spIdx][n][a];
			}
			kineticSumR2 += uR[i] * sD2[spIdx][n];
			kineticSumI2 += uI[i] * sD2[spIdx][n];
			spIdx++;
		}
		int idx = 0;
		for (int i = np2; i < np3; i++)
		{
			for (int a = 0; a < DIM; a++)
			{
				tmp = (bcFactorsEnd[idx][0] * sD[numberOfSplines - 3][n][a] + bcFactorsEnd[idx][1] * sD[numberOfSplines - 2][n][a] + bcFactorsEnd[idx][2] * sD[numberOfSplines - 1][n][a]);
				vecKineticSumR1[a] += uR[i] * tmp;
				vecKineticSumI1[a] += uI[i] * tmp;
			}
			tmp = (bcFactorsEnd[idx][0] * sD2[numberOfSplines - 3][n] + bcFactorsEnd[idx][1] * sD2[numberOfSplines - 2][n] + bcFactorsEnd[idx][2] * sD2[numberOfSplines - 1][n]);
			kineticSumR2 += uR[i] * tmp;
			kineticSumI2 += uI[i] * tmp;
			idx++;
		}

		kineticSumR1I1 += 2.0 * VectorDotProduct_DIM(vecKineticSumR1, vecKineticSumI1);
		kineticSumR1 += VectorNorm2_DIM(vecKineticSumR1);
		kineticSumI1 += VectorNorm2_DIM(vecKineticSumI1);
	}

	kineticR = -(kineticSumR1 - kineticSumI1 + kineticSumR2) * HBAR2_2M;
	kineticI = -(kineticSumR1I1 + kineticSumI2) * HBAR2_2M;

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
	otherExpectationValues[4] = kineticSumR1;
	otherExpectationValues[5] = kineticSumI1;
	otherExpectationValues[6] = kineticSumR2;
	otherExpectationValues[7] = kineticSumI2;
	otherExpectationValues[8] = kineticSumR1I1;
}

void Bosons1D::CalculateExpectationValues(vector<vector<double> >& R, vector<double>& uR, vector<double>& uI, double phiR, double phiI)
{
	RefreshLocalOperators();
	CalculateOtherLocalOperators(R);
	vector<double> dummy_grBins;
	CalculateExpectationValues(this->localOperators, this->splineSumsD, this->splineSumsD2, this->otherLocalOperators, dummy_grBins, uR, uI, phiR, phiI);
}

void Bosons1D::CalculateExpectationValues(ICorrelatedSamplingData* sample, vector<double>& uR, vector<double>& uI, double phiR, double phiI)
{
	CSDataBulkSplines* s = dynamic_cast<CSDataBulkSplines*>(sample);
	CalculateExpectationValues(s->localOperators, s->splineSumsD, s->splineSumsD2, s->otherLocalOperators, s->grBins, uR, uI, phiR, phiI);
}

void Bosons1D::CalculateAdditionalSystemProperties(vector<vector<double> >& R, vector<double>& uR, vector<double>& uI, double phiR, double phiI)
{
	additionalObservables.ClearValues();
	double rni;
	vector<double> vecrni(DIM);

	//pairDistribution
	double weight = 1.0 / ((double)(N - 1));
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < i; j++)
		{
			rni = VectorDisplacementNIC_DIM(R[i], R[j], vecrni);
			if (rni < pairDistribution.grid.max)
			{
				pairDistribution.AddToHistogram(0, rni, weight);
			}
		}
	}

	//structureFactor
	//according to Zhang, Kai. "On the concept of static structure factor." arXiv preprint arXiv:1606.03610 (2016).
	double sumS;
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
				sumS = VectorDotProduct_DIM(kValues[k][kn], R[i]);
				sumSCos[k] += cos(sumS);
				sumSSin[k] += sin(sumS);
			}
		}
	}
	for (int k = 0; k < numOfkValues; k++)
	{
		sk[k] = (sumSCos[k] * sumSCos[k] + sumSSin[k] * sumSSin[k]) / ((double) (N * kValues[k].size()));
		structureFactor.SetValueAtGridIndex(0, k, sk[k]);
	}
}

void Bosons1D::CalculateWavefunction(vector<double>& O, vector<double>& uR, vector<double>& uI, double phiR, double phiI)
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

void Bosons1D::CalculateWavefunction(vector<vector<double> >& R, vector<double>& uR, vector<double>& uI, double phiR, double phiI)
{
	CalculateLocalOperators(R);
	CalculateWavefunction(this->localOperators, uR, uI, phiR, phiI);
}

void Bosons1D::CalculateWavefunction(ICorrelatedSamplingData* sample, vector<double>& uR, vector<double>& uI, double phiR, double phiI)
{
	CSDataBulkSplines* s = dynamic_cast<CSDataBulkSplines*>(sample);
	CalculateWavefunction(s->localOperators, uR, uI, phiR, phiI);
}

void Bosons1D::CalculateWFChange(vector<vector<double> >& R, vector<double>& uR, vector<double>& uI, double phiR, double phiI, int changedParticleIndex, vector<double>& oldPosition)
{
	double sum = 0;
	int bin;
	double rni;
	double rni2;
	double rni3;
	vector<double> vecrni(DIM);

	ClearVector(sumOldPerBin);
	ClearVector(sumNewPerBin);
	for (int i = particleStartIndexForReducedSampling; i < N; i++)
	{
		if (i != changedParticleIndex)
		{
			rni = VectorDisplacementNIC_DIM(R[i], oldPosition, vecrni);
			if (rni <= maxDistance)
			{
				auto interval = lower_bound(nodes.begin(), nodes.end(), rni);
				bin = interval - nodes.begin() - 1;
				if (bin >= numberOfSplines)
				{
					cout << "!!!!" << endl;
				}
				rni2 = rni * rni;
				rni3 = rni2 * rni;

				if (usePreCalcSplineValues)
				{
					int rniBin = (rni - nodes[bin]) * binSizePerNode_R[bin + 3];
					for (int p = 0; p < 4; p++)
					{
						sumOldPerBin[bin - p] += preCalcSplineValues[bin - p][p][rniBin];
					}
				}
				else
				{
					for (int p = 0; p < 4; p++)
					{
						sumOldPerBin[bin - p] += splineWeights[bin - p][p][0] + splineWeights[bin - p][p][1] * rni + splineWeights[bin - p][p][2] * rni2 + splineWeights[bin - p][p][3] * rni3;
					}
				}
			}
			else
			{
				cout << "!!!" << endl;
			}
			rni = VectorDisplacementNIC_DIM(R[i], R[changedParticleIndex], vecrni);
			if (rni <= maxDistance)
			{
				auto interval = lower_bound(nodes.begin(), nodes.end(), rni);
				bin = interval - nodes.begin() - 1;
				rni2 = rni * rni;
				rni3 = rni2 * rni;

				if (usePreCalcSplineValues)
				{
					int rniBin = (rni - nodes[bin]) * binSizePerNode_R[bin + 3];
					for (int p = 0; p < 4; p++)
					{
						sumNewPerBin[bin - p] += preCalcSplineValues[bin - p][p][rniBin];
					}
				}
				else
				{
					for (int p = 0; p < 4; p++)
					{
						sumNewPerBin[bin - p] += splineWeights[bin - p][p][0] + splineWeights[bin - p][p][1] * rni + splineWeights[bin - p][p][2] * rni2 + splineWeights[bin - p][p][3] * rni3;
					}
				}
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
	int spIdx = 3;
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

double Bosons1D::CalculateWFQuotient(vector<vector<double> >& R, vector<double>& uR, vector<double>& uI, double phiR, double phiI, int changedParticleIndex, vector<double>& oldPosition)
{
	CalculateWFChange(R, uR, uI, phiR, phiI, changedParticleIndex, oldPosition);
	double wfQuotient = exp(2.0 * (exponentNew - exponent));

	//cout << "qu=" << wfQuotient << " (" << exponentNew << " - " << exponent << ")" << endl;

	return wfQuotient;
}

void Bosons1D::AcceptMove()
{
	wf = wfNew;
	exponent = exponentNew;
	splineSums = splineSumsNew;
}

void Bosons1D::InitCorrelatedSamplingData(vector<ICorrelatedSamplingData*>& data)
{
	for (unsigned int i = 0; i < data.size(); i++)
	{
		data[i] = new CSDataBulkSplines();
	}
}

void Bosons1D::FillCorrelatedSamplingData(ICorrelatedSamplingData* data)
{
	CSDataBulkSplines* d = dynamic_cast<CSDataBulkSplines*>(data);
	d->splineSumsD = this->splineSumsD;
	d->splineSumsD2 = this->splineSumsD2;
	d->otherLocalOperators = this->otherLocalOperators;
}

}
