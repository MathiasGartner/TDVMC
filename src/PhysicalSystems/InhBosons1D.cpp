#include "InhBosons1D.h"

#include "../SplineFactory.h"

namespace PhysicalSystems
{

InhBosons1D::InhBosons1D(vector<double>& params, string configDirectory) :
		IPhysicalSystem(params, configDirectory)
{
	this->USE_NORMALIZATION_AND_PHASE = true;
	this->USE_NIC = true;
	this->USE_MOVE_COM_TO_ZERO = false;

	halfLength = 0;
	maxDistance = 0;

	numOfOtherLocalOperators = 0;
	numOfkValues = 0;

	exponentNew = 0;
}

void InhBosons1D::ExtendNodes(vector<double>& n)
{
	int periodicNodeCount = 3;
	for (int i = 0; i < periodicNodeCount; i++)
	{
		n.insert(n.begin(), -n[2 * i + 1]);
		n.push_back(2.0 * n[n.size() - 1 - i] - n[n.size() - 2 * (i + 1)]);
	}
}

void InhBosons1D::SetNodesSPF(vector<double> n)
{
	this->spf.nodes = n;
	ExtendNodes(this->spf.nodes);
}

void InhBosons1D::SetNodesPC(vector<double> n)
{
	this->pc.nodes = n;
	ExtendNodes(this->pc.nodes);
}

void InhBosons1D::InitSystem()
{
	numOfOtherLocalOperators = 4; //INFO: potentialIntern, potentialInternComplex, potentialExtern, potentialExternComplex
	otherLocalOperators.resize(numOfOtherLocalOperators);
	numOfkValues = 50;
	numOfOtherExpectationValues = 9;
	numOfAdditionalSystemProperties = numOfOtherExpectationValues;

	halfLength = LBOX / 2.0;
	maxDistance = pc.nodes[pc.nodes.size() - 4]; //INFO: only valid for u_2(rij) that is defined up to L/2

	//INFO: use same nodes and spf defined in the interval [0, L/2] only for symmetric external potentials
	spf.numberOfSplines = spf.nodes.size() - (3 + 1); //3rd order splines -(d+1)
	spf.splineWeights = SplineFactory::GetWeights(spf.nodes);
	SplineFactory::SetBoundaryConditions3_1D_OR_2(spf.nodes, spf.bcFactorsStart, true);
	SplineFactory::SetBoundaryConditions3_1D_CO_2(spf.nodes, spf.bcFactorsEnd, maxDistance, true);
	spf.numberOfSpecialParametersStart = spf.bcFactorsStart.size();
	spf.numberOfSpecialParametersEnd = spf.bcFactorsEnd.size();
	spf.numberOfStandardParameters = spf.nodes.size() - 2 * spf.splineOrder - spf.numberOfSpecialParametersStart - spf.numberOfSpecialParametersEnd;
	spf.np1 = spf.numberOfSpecialParametersStart;
	spf.np2 = spf.np1 + spf.numberOfStandardParameters;
	spf.np3 = spf.np2 + spf.numberOfSpecialParametersEnd;

	pc.numberOfSplines = pc.nodes.size() - (3 + 1); //3rd order splines -(d+1)
	pc.splineWeights = SplineFactory::GetWeights(pc.nodes);
	SplineFactory::SetBoundaryConditions3_1D_OR_2(pc.nodes, pc.bcFactorsStart, true);
	SplineFactory::SetBoundaryConditions3_1D_CO_2(pc.nodes, pc.bcFactorsEnd, maxDistance, true);
	pc.numberOfSpecialParametersStart = pc.bcFactorsStart.size();
	pc.numberOfSpecialParametersEnd = pc.bcFactorsEnd.size();
	pc.numberOfStandardParameters = pc.nodes.size() - 2 * pc.splineOrder - pc.numberOfSpecialParametersStart - pc.numberOfSpecialParametersEnd;
	pc.np1 = pc.numberOfSpecialParametersStart;
	pc.np2 = pc.np1 + pc.numberOfStandardParameters;
	pc.np3 = pc.np2 + pc.numberOfSpecialParametersEnd;

	int requiredParams = spf.numberOfSplines - (3 - spf.numberOfSpecialParametersStart) - (3 - spf.numberOfSpecialParametersEnd) +
						 pc.numberOfSplines - (3 - pc.numberOfSpecialParametersStart) - (3 - pc.numberOfSpecialParametersEnd);
	if (N_PARAM != requiredParams)
	{
		cout << "!!! WRONG NUMBER OF PARAMETERS !!!" << endl;
		cout << "!!! need " << requiredParams << " parameters, got " << N_PARAM << " !!!" << endl;
		exit(0);
	}

	wf = 0.0;
	exponent = 0.0;
	InitVector(spf.splineSums, spf.numberOfSplines, 0.0);
	InitVector(spf.splineSumsD, spf.numberOfSplines, N, DIM, 0.0);
	InitVector(spf.splineSumsD2, spf.numberOfSplines, N, 0.0);
	InitVector(pc.splineSums, pc.numberOfSplines, 0.0);
	InitVector(pc.splineSumsD, pc.numberOfSplines, N, DIM, 0.0);
	InitVector(pc.splineSumsD2, pc.numberOfSplines, N, 0.0);

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
	InitVector(spf.splineSumsNew, spf.numberOfSplines, 0.0);
	InitVector(spf.sumOldPerBin, spf.numberOfSplines, 0.0);
	InitVector(spf.sumNewPerBin, spf.numberOfSplines, 0.0);
	InitVector(pc.splineSumsNew, pc.numberOfSplines, 0.0);
	InitVector(pc.sumOldPerBin, pc.numberOfSplines, 0.0);
	InitVector(pc.sumNewPerBin, pc.numberOfSplines, 0.0);

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

	density.name = "density";
	density.grid.name = "r";
	density.InitGrid(0.0, halfLength, halfLength / 100.0);
	density.InitScaling();
	density.InitObservables( { "rho(r)" });

	pairDistribution.name = "pairDistribution";
	pairDistribution.grid.name = "r_ij";
	//pairDistribution.InitGrid(0.0, halfLength, 0.05);
	//pairDistribution.InitGrid(0.0, pc.nodes[pc.nodes.size() - 4], pc.nodes[pc.nodes.size() - 4] / 100.0);
	pairDistribution.InitGrid(0.0, halfLength, halfLength / 100.0);
	pairDistribution.InitScaling();
	pairDistribution.InitObservables( { "g_2(r_ij)" });

	structureFactor.name = "structureFactor";
	structureFactor.grid.name = "k";
	structureFactor.InitGrid(kNorms);
	structureFactor.InitObservables( { "S(k)" });

	additionalObservables.Add(&density);
	additionalObservables.Add(&pairDistribution);
	additionalObservables.Add(&structureFactor);
}

void InhBosons1D::RefreshLocalOperators()
{
	int offset;
	int spIdx;
	int idx;

	//INFO: SingleParticleFunction
	offset = 0;
	for (int i = 0; i < spf.np1; i++)
	{
		this->localOperators[i] = (spf.bcFactorsStart[i][0] * spf.splineSums[0] + spf.bcFactorsStart[i][1] * spf.splineSums[1] + spf.bcFactorsStart[i][2] * spf.splineSums[2]);
	}
	spIdx = 3;
	for (int i = spf.np1; i < spf.np2; i++)
	{
		this->localOperators[i] = spf.splineSums[spIdx];
		spIdx++;
	}
	idx = 0;
	for (int i = spf.np2; i < spf.np3; i++)
	{
		this->localOperators[i] = (spf.bcFactorsEnd[idx][0] * spf.splineSums[spf.numberOfSplines - 3] + spf.bcFactorsEnd[idx][1] * spf.splineSums[spf.numberOfSplines - 2] + spf.bcFactorsEnd[idx][2] * spf.splineSums[spf.numberOfSplines - 1]);
		idx++;
	}

	//INFO: PairCorrelation
	offset = spf.np3;
	for (int i = 0; i < pc.np1; i++)
	{
		this->localOperators[offset + i] = (pc.bcFactorsStart[i][0] * pc.splineSums[0] + pc.bcFactorsStart[i][1] * pc.splineSums[1] + pc.bcFactorsStart[i][2] * pc.splineSums[2]);
	}
	spIdx = 3;
	for (int i = pc.np1; i < pc.np2; i++)
	{
		this->localOperators[offset + i] = pc.splineSums[spIdx];
		spIdx++;
	}
	idx = 0;
	for (int i = pc.np2; i < pc.np3; i++)
	{
		this->localOperators[offset + i] = (pc.bcFactorsEnd[idx][0] * pc.splineSums[pc.numberOfSplines - 3] + pc.bcFactorsEnd[idx][1] * pc.splineSums[pc.numberOfSplines - 2] + pc.bcFactorsEnd[idx][2] * pc.splineSums[pc.numberOfSplines - 1]);
		idx++;
	}
}

void InhBosons1D::CalculateLocalOperators(vector<vector<double> >& R)
{
	int bin;
	double rni;
	double rni2;
	double rni3;
	double r;
	double r2;
	double r3;
	vector<double> vecrni(DIM);

	ClearVector(spf.splineSums);
	ClearVector(pc.splineSums);

	for (int n = 0; n < N; n++)
	{
		//INFO: SingleParticleFunction
		r = abs(GetCoordinateNIC(R[n][0])); //INFO: only for 1D!
		if (r <= maxDistance)
		{
			bin = GetBinIndex(spf.nodes, r);
			r2 = r * r;
			r3 = r2 * r;
			for (int p = 0; p < spf.numOfSplineParts; p++)
			{
				spf.splineSums[bin - p] += spf.splineWeights[bin - p][p][0] + spf.splineWeights[bin - p][p][1] * r + spf.splineWeights[bin - p][p][2] * r2 + spf.splineWeights[bin - p][p][3] * r3;
			}
		}

		//INFO: PairCorrelation
		for (int i = 0; i < n; i++)
		{
			rni = VectorDisplacementNIC_DIM(R[n], R[i], vecrni);
			if (rni <= maxDistance)
			{
				bin = GetBinIndex(pc.nodes, rni);
				rni2 = rni * rni;
				rni3 = rni2 * rni;

				for (int p = 0; p < pc.numOfSplineParts; p++)
				{
					pc.splineSums[bin - p] += pc.splineWeights[bin - p][p][0] + pc.splineWeights[bin - p][p][1] * rni + pc.splineWeights[bin - p][p][2] * rni2 + pc.splineWeights[bin - p][p][3] * rni3;
				}
			}
			else
			{
				cout << "!!!!" << endl;
			}
		}
	}
	RefreshLocalOperators();
}

void InhBosons1D::CalculateOtherLocalOperators(vector<vector<double> >& R)
{
	vector<double> vecr(DIM);
	vector<double> evecr(DIM);
	vector<double> vecrni(DIM);
	vector<double> evecrni(DIM);
	vector<double> spftmp1(spf.numOfSplineParts);
	vector<double> spftmp2(spf.numOfSplineParts);
	vector<double> pctmp1(pc.numOfSplineParts);
	vector<double> pctmp2(pc.numOfSplineParts);
	int bin;
	double rni;
	double rni2;
	double r;
	double r2;

	double potentialExtern = 0;
	double potentialExternComplex = 0;
	double potentialIntern = 0;
	double potentialInternComplex = 0;

	//potential parameters (a -> width; b -> height)
	double potentialRange = params[0];
	double potentialStrength = params[1];
	if (params.size() > 4 && this->time >= params[4])
	{
		potentialRange = params[5];
		potentialStrength = params[6];
	}

	ClearVector(spf.splineSumsD);
	ClearVector(spf.splineSumsD2);
	ClearVector(pc.splineSumsD);
	ClearVector(pc.splineSumsD2);
	ClearVector(otherLocalOperators);

	for (int n = 0; n < N; n++)
	{
		//external potential energy
		potentialExtern += GetExternalPotential(R[n]);

		//INFO: SingleParticleFunction
		r = abs(GetCoordinateNIC(R[n][0])); //INFO: only for 1D!
		bin = GetBinIndex(spf.nodes, r);
		r2 = r * r;

		for (int p = 0; p < spf.numOfSplineParts; p++)
		{
			//first derivative
			spftmp1[spf.numOfSplineParts - 1 - p] = spf.splineWeights[bin - p][p][1] + 2.0 * spf.splineWeights[bin - p][p][2] * r + 3.0 * spf.splineWeights[bin - p][p][3] * r2;
			//second derivative
			spftmp2[spf.numOfSplineParts - 1 - p] = 2.0 * spf.splineWeights[bin - p][p][2] + 6.0 * spf.splineWeights[bin - p][p][3] * r;
		}


		evecr[0] = GetCoordinateNIC(R[n][0]) > 0 ? 1.0 : -1.0; //TODO: is this the unit vector? do not calculate again
		bool outsideL2 = false;
		double firstDerivativeFactor = outsideL2 ? -1.0 : 1.0;
		for (int a = 0; a < DIM; a++)
		{
			for (int b = 0; b < spf.numOfSplineParts; b++)
			{
				spf.splineSumsD[bin - b][n][a] += spftmp1[spf.numOfSplineParts - 1 - b] * evecr[a] * firstDerivativeFactor;
			}
		}
		double secondDerivativeFactor = DIM - 1.0; //INFO: at least correct for 2D and 3D (and 1D?)
		for (int b = 0; b < spf.numOfSplineParts; b++)
		{
			spf.splineSumsD2[bin - b][n] += spftmp2[spf.numOfSplineParts - 1 - b] + secondDerivativeFactor / r * spftmp1[spf.numOfSplineParts - 1 - b];
		}

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
					if (rni < potentialRange)
					{
						potentialIntern += potentialStrength;
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
					bin = GetBinIndex(pc.nodes, rni);
					rni2 = rni * rni;

					for (int p = 0; p < pc.numOfSplineParts; p++)
					{
						//first derivative
						pctmp1[pc.numOfSplineParts - 1 - p] = pc.splineWeights[bin - p][p][1] + 2.0 * pc.splineWeights[bin - p][p][2] * rni + 3.0 * pc.splineWeights[bin - p][p][3] * rni2;
						//second derivative
						pctmp2[pc.numOfSplineParts - 1 - p] = 2.0 * pc.splineWeights[bin - p][p][2] + 6.0 * pc.splineWeights[bin - p][p][3] * rni;
					}

					for (int a = 0; a < DIM; a++)
					{
						evecrni[a] = vecrni[a] / rni;
					}
					bool outsideL2 = false;
					double firstDerivativeFactor = outsideL2 ? -1.0 : 1.0;
					for (int a = 0; a < DIM; a++)
					{
						for (int b = 0; b < pc.numOfSplineParts; b++)
						{
							pc.splineSumsD[bin - b][n][a] += pctmp1[pc.numOfSplineParts - 1 - b] * evecrni[a] * firstDerivativeFactor;
						}
					}
					double secondDerivativeFactor = DIM - 1.0; //INFO: at least correct for 2D and 3D (and 1D?)
					for (int b = 0; b < pc.numOfSplineParts; b++)
					{
						pc.splineSumsD2[bin - b][n] += pctmp2[pc.numOfSplineParts - 1 - b] + secondDerivativeFactor / rni * pctmp1[pc.numOfSplineParts - 1 - b];
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

vector<double> InhBosons1D::GetCenterOfMass(vector<vector<double> >& R)
{
	vector<double> com = { 0.0, 0.0, 0.0 };
	return com;
}

double InhBosons1D::GetExternalPotential(vector<double>& r)
{
	double value;
	double potentialFrequency = params[2];
	double potentialStrength = params[3];
	if (params.size() > 4 && this->time >= params[4])
	{
		potentialFrequency = params[7];
		potentialStrength = params[8];
	}
	value = cos(2.0 * M_PI * r[0] / (maxDistance / potentialFrequency));
	value *= potentialStrength;
	//value = 0.0;
	return value;
}

void InhBosons1D::CalculateExpectationValues(vector<double>& O, vector<vector<vector<double> > >& spf_sD, vector<vector<double> >& spf_sD2, vector<vector<vector<double> > >& pc_sD, vector<vector<double> >& pc_sD2, vector<double>& otherO, vector<double>& uR, vector<double>& uI, double phiR, double phiI)
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
	int offset;
	int spIdx;
	int idx;

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

		//INFO: SingleParticleFunction
		offset = 0;
		for (int i = 0; i < spf.np1; i++)
		{
			for (int a = 0; a < DIM; a++)
			{
				tmp = (spf.bcFactorsStart[i][0] * spf_sD[0][n][a] + spf.bcFactorsStart[i][1] * spf_sD[1][n][a] + spf.bcFactorsStart[i][2] * spf_sD[2][n][a]);
				vecKineticSumR1[a] += uR[i] * tmp;
				vecKineticSumI1[a] += uI[i] * tmp;
			}
			tmp = (spf.bcFactorsStart[i][0] * spf_sD2[0][n] + spf.bcFactorsStart[i][1] * spf_sD2[1][n] + spf.bcFactorsStart[i][2] * spf_sD2[2][n]);
			kineticSumR2 += uR[i] * tmp;
			kineticSumI2 += uI[i] * tmp;
		}
		spIdx = 3;
		for (int i = spf.np1; i < spf.np2; i++)
		{
			for (int a = 0; a < DIM; a++)
			{
				vecKineticSumR1[a] += uR[i] * spf_sD[spIdx][n][a];
				vecKineticSumI1[a] += uI[i] * spf_sD[spIdx][n][a];
			}
			kineticSumR2 += uR[i] * spf_sD2[spIdx][n];
			kineticSumI2 += uI[i] * spf_sD2[spIdx][n];
			spIdx++;
		}
		idx = 0;
		for (int i = spf.np2; i < spf.np3; i++)
		{
			for (int a = 0; a < DIM; a++)
			{
				tmp = (spf.bcFactorsEnd[idx][0] * spf_sD[spf.numberOfSplines - 3][n][a] + spf.bcFactorsEnd[idx][1] * spf_sD[spf.numberOfSplines - 2][n][a] + spf.bcFactorsEnd[idx][2] * spf_sD[spf.numberOfSplines - 1][n][a]);
				vecKineticSumR1[a] += uR[i] * tmp;
				vecKineticSumI1[a] += uI[i] * tmp;
			}
			tmp = (spf.bcFactorsEnd[idx][0] * spf_sD2[spf.numberOfSplines - 3][n] + spf.bcFactorsEnd[idx][1] * spf_sD2[spf.numberOfSplines - 2][n] + spf.bcFactorsEnd[idx][2] * spf_sD2[spf.numberOfSplines - 1][n]);
			kineticSumR2 += uR[i] * tmp;
			kineticSumI2 += uI[i] * tmp;
			idx++;
		}

		//INFO: PairCorrelation
		offset = spf.np3;
		for (int i = 0; i < pc.np1; i++)
		{
			for (int a = 0; a < DIM; a++)
			{
				tmp = (pc.bcFactorsStart[i][0] * pc_sD[0][n][a] + pc.bcFactorsStart[i][1] * pc_sD[1][n][a] + pc.bcFactorsStart[i][2] * pc_sD[2][n][a]);
				vecKineticSumR1[a] += uR[offset + i] * tmp;
				vecKineticSumI1[a] += uI[offset + i] * tmp;
			}
			tmp = (pc.bcFactorsStart[i][0] * pc_sD2[0][n] + pc.bcFactorsStart[i][1] * pc_sD2[1][n] + pc.bcFactorsStart[i][2] * pc_sD2[2][n]);
			kineticSumR2 += uR[offset + i] * tmp;
			kineticSumI2 += uI[offset + i] * tmp;
		}
		spIdx = 3;
		for (int i = pc.np1; i < pc.np2; i++)
		{
			for (int a = 0; a < DIM; a++)
			{
				vecKineticSumR1[a] += uR[offset + i] * pc_sD[spIdx][n][a];
				vecKineticSumI1[a] += uI[offset + i] * pc_sD[spIdx][n][a];
			}
			kineticSumR2 += uR[offset + i] * pc_sD2[spIdx][n];
			kineticSumI2 += uI[offset + i] * pc_sD2[spIdx][n];
			spIdx++;
		}
		idx = 0;
		for (int i = pc.np2; i < pc.np3; i++)
		{
			for (int a = 0; a < DIM; a++)
			{
				tmp = (pc.bcFactorsEnd[idx][0] * pc_sD[pc.numberOfSplines - 3][n][a] + pc.bcFactorsEnd[idx][1] * pc_sD[pc.numberOfSplines - 2][n][a] + pc.bcFactorsEnd[idx][2] * pc_sD[pc.numberOfSplines - 1][n][a]);
				vecKineticSumR1[a] += uR[offset + i] * tmp;
				vecKineticSumI1[a] += uI[offset + i] * tmp;
			}
			tmp = (pc.bcFactorsEnd[idx][0] * pc_sD2[pc.numberOfSplines - 3][n] + pc.bcFactorsEnd[idx][1] * pc_sD2[pc.numberOfSplines - 2][n] + pc.bcFactorsEnd[idx][2] * pc_sD2[pc.numberOfSplines - 1][n]);
			kineticSumR2 += uR[offset + i] * tmp;
			kineticSumI2 += uI[offset + i] * tmp;
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

void InhBosons1D::CalculateExpectationValues(vector<vector<double> >& R, vector<double>& uR, vector<double>& uI, double phiR, double phiI)
{
	RefreshLocalOperators();
	CalculateOtherLocalOperators(R);
	CalculateExpectationValues(this->localOperators, this->spf.splineSumsD, this->spf.splineSumsD2, this->pc.splineSumsD, this->pc.splineSumsD2, this->otherLocalOperators, uR, uI, phiR, phiI);
}

void InhBosons1D::CalculateExpectationValues(ICorrelatedSamplingData* sample, vector<double>& uR, vector<double>& uI, double phiR, double phiI)
{
	cout << "not implemented" << endl;
}

void InhBosons1D::CalculateAdditionalSystemProperties(vector<vector<double> >& R, vector<double>& uR, vector<double>& uI, double phiR, double phiI)
{
	additionalObservables.ClearValues();
	double r;
	double rni;
	vector<double> vecrni(DIM);

	//density
	for (int i = 0; i < N; i++)
	{
		r = abs(GetCoordinateNIC(R[i][0]));
		if (r < density.grid.max)
		{
			density.AddToHistogram(0, r, 1.0); //INFO: only for 1D - otherwise there should be a multi-dim histogram variable r
		}
	}

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

void InhBosons1D::CalculateWavefunction(vector<double>& O, vector<double>& uR, vector<double>& uI, double phiR, double phiI)
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

void InhBosons1D::CalculateWavefunction(vector<vector<double> >& R, vector<double>& uR, vector<double>& uI, double phiR, double phiI)
{
	CalculateLocalOperators(R);
	CalculateWavefunction(this->localOperators, uR, uI, phiR, phiI);
}

void InhBosons1D::CalculateWavefunction(ICorrelatedSamplingData* sample, vector<double>& uR, vector<double>& uI, double phiR, double phiI)
{
	CSDataBulkSplines* s = dynamic_cast<CSDataBulkSplines*>(sample);
	CalculateWavefunction(s->localOperators, uR, uI, phiR, phiI);
}

void InhBosons1D::CalculateWFChange(vector<vector<double> >& R, vector<double>& uR, vector<double>& uI, double phiR, double phiI, int changedParticleIndex, vector<double>& oldPosition)
{
	double sum = 0;
	int bin = 0;
	double rni;
	double rni2;
	double rni3;
	double r;
	double r2;
	double r3;
	vector<double> vecrni(DIM);

	ClearVector(spf.sumOldPerBin);
	ClearVector(spf.sumNewPerBin);
	ClearVector(pc.sumOldPerBin);
	ClearVector(pc.sumNewPerBin);

	//INFO: SingleParticleFunction
	r = abs(GetCoordinateNIC(oldPosition[0])); //INFO: only for 1D!
	if (r <= maxDistance)
	{
		bin = GetBinIndex(spf.nodes, r);
		r2 = r * r;
		r3 = r2 * r;
		for (int p = 0; p < spf.numOfSplineParts; p++)
		{
			spf.sumOldPerBin[bin - p] += spf.splineWeights[bin - p][p][0] + spf.splineWeights[bin - p][p][1] * r + spf.splineWeights[bin - p][p][2] * r2 + spf.splineWeights[bin - p][p][3] * r3;
		}
	}
	r = abs(GetCoordinateNIC(R[changedParticleIndex][0])); //INFO: only for 1D!
	if (r <= maxDistance)
	{
		bin = GetBinIndex(spf.nodes, r);
		r2 = r * r;
		r3 = r2 * r;
		for (int p = 0; p < spf.numOfSplineParts; p++)
		{
			spf.sumNewPerBin[bin - p] += spf.splineWeights[bin - p][p][0] + spf.splineWeights[bin - p][p][1] * r + spf.splineWeights[bin - p][p][2] * r2 + spf.splineWeights[bin - p][p][3] * r3;
		}
	}

	//INFO: PairCorrelation
	for (int i = 0; i < N; i++)
	{
		if (i != changedParticleIndex)
		{
			rni = VectorDisplacementNIC_DIM(R[i], oldPosition, vecrni);
			if (rni <= maxDistance)
			{
				bin = GetBinIndex(pc.nodes, rni);
				rni2 = rni * rni;
				rni3 = rni2 * rni;

				for (int p = 0; p < pc.numOfSplineParts; p++)
				{
					pc.sumOldPerBin[bin - p] += pc.splineWeights[bin - p][p][0] + pc.splineWeights[bin - p][p][1] * rni + pc.splineWeights[bin - p][p][2] * rni2 + pc.splineWeights[bin - p][p][3] * rni3;
				}
			}
			else
			{
				cout << "!!!" << endl;
			}
			rni = VectorDisplacementNIC_DIM(R[i], R[changedParticleIndex], vecrni);
			if (rni <= maxDistance)
			{
				bin = GetBinIndex(pc.nodes, rni);
				rni2 = rni * rni;
				rni3 = rni2 * rni;

				for (int p = 0; p < pc.numOfSplineParts; p++)
				{
					pc.sumNewPerBin[bin - p] += pc.splineWeights[bin - p][p][0] + pc.splineWeights[bin - p][p][1] * rni + pc.splineWeights[bin - p][p][2] * rni2 + pc.splineWeights[bin - p][p][3] * rni3;
				}
			}
			else
			{
				cout << "!!!" << endl;
			}
		}
	}

	//INFO: check that there is no negative sum due to rounding/floating point errors
	for (int i = 0; i < spf.numberOfSplines; i++)
	{
		spf.splineSumsNew[i] = max(0.0, spf.splineSums[i] - spf.sumOldPerBin[i] + spf.sumNewPerBin[i]);
	}
	for (int i = 0; i < pc.numberOfSplines; i++)
	{
		pc.splineSumsNew[i] = max(0.0, pc.splineSums[i] - pc.sumOldPerBin[i] + pc.sumNewPerBin[i]);
	}


	int offset;
	int spIdx;
	int idx;

	//INFO: SingleParticleFunction
	offset = 0;
	for (int i = 0; i < spf.np1; i++)
	{
		sum += uR[i] * (spf.bcFactorsStart[i][0] * spf.splineSumsNew[0] + spf.bcFactorsStart[i][1] * spf.splineSumsNew[1] + spf.bcFactorsStart[i][2] * spf.splineSumsNew[2]);
	}
	spIdx = 3;
	for (int i = spf.np1; i < spf.np2; i++)
	{
		sum += uR[i] * spf.splineSumsNew[spIdx];
		spIdx++;
	}
	idx = 0;
	for (int i = spf.np2; i < spf.np3; i++)
	{
		sum += uR[i] * (spf.bcFactorsEnd[idx][0] * spf.splineSumsNew[spf.numberOfSplines - 3] + spf.bcFactorsEnd[idx][1] * spf.splineSumsNew[spf.numberOfSplines - 2] + spf.bcFactorsEnd[idx][2] * spf.splineSumsNew[spf.numberOfSplines - 1]);
		idx++;
	}

	//INFO: PairCorrelation
	offset = spf.np3;
	for (int i = 0; i < pc.np1; i++)
	{
		sum += uR[offset + i] * (pc.bcFactorsStart[i][0] * pc.splineSumsNew[0] + pc.bcFactorsStart[i][1] * pc.splineSumsNew[1] + pc.bcFactorsStart[i][2] * pc.splineSumsNew[2]);
	}
	spIdx = 3;
	for (int i = pc.np1; i < pc.np2; i++)
	{
		sum += uR[offset + i] * pc.splineSumsNew[spIdx];
		spIdx++;
	}
	idx = 0;
	for (int i = pc.np2; i < pc.np3; i++)
	{
		sum += uR[offset + i] * (pc.bcFactorsEnd[idx][0] * pc.splineSumsNew[pc.numberOfSplines - 3] + pc.bcFactorsEnd[idx][1] * pc.splineSumsNew[pc.numberOfSplines - 2] + pc.bcFactorsEnd[idx][2] * pc.splineSumsNew[pc.numberOfSplines - 1]);
		idx++;
	}

	exponentNew = sum;
	wfNew = exp(exponentNew + phiR);
}

double InhBosons1D::CalculateWFQuotient(vector<vector<double> >& R, vector<double>& uR, vector<double>& uI, double phiR, double phiI, int changedParticleIndex, vector<double>& oldPosition)
{
	CalculateWFChange(R, uR, uI, phiR, phiI, changedParticleIndex, oldPosition);
	double wfQuotient = exp(2.0 * (exponentNew - exponent));

	//cout << "qu=" << wfQuotient << " (" << exponentNew << " - " << exponent << ")" << endl;

	return wfQuotient;
}

void InhBosons1D::AcceptMove()
{
	wf = wfNew;
	exponent = exponentNew;
	spf.splineSums = spf.splineSumsNew;
	pc.splineSums = pc.splineSumsNew;
}

void InhBosons1D::InitCorrelatedSamplingData(vector<ICorrelatedSamplingData*>& data)
{
}

void InhBosons1D::FillCorrelatedSamplingData(ICorrelatedSamplingData* data)
{
}

}
