#include "InhContactBosons.h"

#include "../SplineFactory.h"

namespace PhysicalSystems
{

InhContactBosons::InhContactBosons(vector<double>& params, string configDirectory) :
		IPhysicalSystem(params, configDirectory)
{
	this->USE_NORMALIZATION_AND_PHASE = true;
	this->USE_NIC = true;
	this->USE_MOVE_COM_TO_ZERO = false;

	halfLength = 0;
	maxDistance = 0;

	numOfOtherLocalOperators = 0;
	numOfPairDistributionValues = 0;
	numOfDensityValues = 0;
	numOfkValues = 0;

	exponentNew = 0;

	//INFO: gamma=10 equals c=20 in 10.1103/PhysRevA.74.031605
	gamma = params[0] == 0.0 ? params[1] : 0.0;
	double potentialK = params[2];
	double kf = M_PI;
	gamma *= potentialK * kf;
}

void InhContactBosons::ExtendNodes(vector<double>& n)
{
	int periodicNodeCount = 3;
	for (int i = 0; i < periodicNodeCount; i++)
	{
		n.insert(n.begin(), -n[2 * i + 1]);
		n.push_back(2.0 * n[n.size() - 1 - i] - n[n.size() - 2 * (i + 1)]);
	}
}

void InhContactBosons::SetNodesSPF(vector<double> n)
{
	this->spf.nodes = n;
	ExtendNodes(this->spf.nodes);
}

void InhContactBosons::SetNodesPC(vector<double> n)
{
	this->pc.nodes = n;
	ExtendNodes(this->pc.nodes);
}

void InhContactBosons::SetPairDistributionBinCount(double n)
{
	this->numOfPairDistributionValues = n;
}

void InhContactBosons::SetDensityBinCount(double n)
{
	this->numOfDensityValues = n;
}

void InhContactBosons::InitSystem()
{
	numOfOtherLocalOperators = 4; //INFO: potentialIntern, potentialInternComplex, potentialExtern, potentialExternComplex
	otherLocalOperators.resize(numOfOtherLocalOperators);
	numOfkValues = 100;
	numOfOtherExpectationValues = 9;
	numOfAdditionalSystemProperties = numOfOtherExpectationValues;

	halfLength = LBOX / 2.0;

	if (spf.nodes.empty())
	{
		//INFO: use half of the parameters for spf
		int nparam_2 = N_PARAM / 2;
		spf.nodeSpacing = LBOX / (double) nparam_2;
		for (int i = -3; i < nparam_2 + 4; i++)
		{
			spf.nodes.push_back((i * LBOX) / ((double) nparam_2));
		}
	}

	if (pc.nodes.empty())
	{
		//INFO: use half of the parameters for pc
		//INFO: BC
		//this is for SetBoundaryConditions3_1D_OR_2 and SetBoundaryConditions3_1D_CO_2
		//and an even number of parameters
		int nparam_2 = N_PARAM / 2;
		pc.nodeSpacing = halfLength / ((double) nparam_2 - 1.0);
		for (int i = -3; i < nparam_2 + 3; i++)
		{
			pc.nodes.push_back((i * halfLength) / ((double) nparam_2 - 1.0));
		}
	}

	maxDistance = pc.nodes[pc.nodes.size() - 4]; //INFO: only valid for u_2(rij) that is defined up to L/2

	spf.numberOfSplines = spf.nodes.size() - (3 + 1); //3rd order splines -(d+1)
	spf.splineWeights = SplineFactory::GetWeights(spf.nodes);
	SplineFactory::SetBoundaryConditions3_1D_OR_1(spf.nodes, spf.bcFactorsStart, true);
	SplineFactory::SetBoundaryConditions3_1D_CO_1(spf.nodes, spf.bcFactorsEnd, LBOX, true);
	spf.numberOfSpecialParametersStart = spf.bcFactorsStart.size();
	spf.numberOfSpecialParametersEnd = spf.bcFactorsEnd.size();
	spf.numberOfStandardParameters = spf.nodes.size() - 2 * spf.splineOrder - spf.numberOfSpecialParametersStart - 1;
	spf.np1 = spf.numberOfSpecialParametersStart;
	spf.np2 = spf.np1 + spf.numberOfStandardParameters;
	spf.np3 = spf.np2; // + spf.numberOfSpecialParametersEnd;

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

	int requiredParams = spf.numberOfSplines - (3 - spf.numberOfSpecialParametersStart) - (3 - spf.numberOfSpecialParametersEnd) - 3 + //INFO: last -3 is because the last three parameters are identical to the first three parameters
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
	density.InitGrid(0.0, LBOX, LBOX / numOfDensityValues);
	density.InitScaling();
	density.InitObservables( { "rho(r)" });

	pairDensity.name = "pairDensity";
	pairDensity.InitGrid( { { 0.0, LBOX, LBOX / numOfPairDistributionValues }, { 0.0, LBOX, LBOX / numOfPairDistributionValues } });
	pairDensity.grids[0].name = "r_i";
	pairDensity.grids[1].name = "r_j";
	pairDensity.InitObservables( { "rho_2(r_i,r_j)" });

	pairDistribution.name = "pairDistribution";
	pairDistribution.grid.name = "r_ij";
	//pairDistribution.InitGrid(0.0, halfLength, 0.05);
	//pairDistribution.InitGrid(0.0, pc.nodes[pc.nodes.size() - 4], pc.nodes[pc.nodes.size() - 4] / 100.0);
	pairDistribution.InitGrid(0.0, halfLength, halfLength / numOfPairDistributionValues);
	pairDistribution.InitScaling();
	pairDistribution.InitObservables( { "g_2(r_ij)" });

	structureFactor.name = "structureFactor";
	structureFactor.grid.name = "k";
	structureFactor.InitGrid(kNorms);
	structureFactor.InitObservables( { "S(k)" });

	additionalObservables.Add(&density);
	additionalObservables.Add(&pairDistribution);
	additionalObservables.Add(&structureFactor);
	additionalObservables.Add(&pairDensity);
}

void InhContactBosons::RefreshLocalOperators()
{
	int offset;
	int spIdx;
	int idx;

	//INFO: SingleParticleFunction
	offset = 0;
	for (int i = 0; i < spf.np1; i++)
	{
		this->localOperators[i] = (spf.bcFactorsStart[i][0] * spf.splineSums[0] + spf.bcFactorsStart[i][1] * spf.splineSums[1] + spf.bcFactorsStart[i][2] * spf.splineSums[2]);
		//INFO: add the last three splineSums to the first three operators
		this->localOperators[i] += (spf.bcFactorsEnd[i][0] * spf.splineSums[spf.numberOfSplines - 3] + spf.bcFactorsEnd[i][1] * spf.splineSums[spf.numberOfSplines - 2] + spf.bcFactorsEnd[i][2] * spf.splineSums[spf.numberOfSplines - 1]);
	}
	spIdx = 3;
	for (int i = spf.np1; i < spf.np2; i++)
	{
		this->localOperators[i] = spf.splineSums[spIdx];
		spIdx++;
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

void InhContactBosons::CalculateLocalOperators(vector<vector<double> >& R)
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
		r = GetCoordinateNIC(R[n][0]) + halfLength; //to get coordinate in interval [0, LBOX]
		if (r < 0 || r > LBOX)
		{
			cout << "!!!!!!";
		}
		else
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

void InhContactBosons::CalculateOtherLocalOperators(vector<vector<double> >& R)
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
		//TODO: set gamma if interaction is quenched
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
		r = GetCoordinateNIC(R[n][0]) + halfLength; //to get coordinate in interval [0, LBOX]
		bin = GetBinIndex(spf.nodes, r);
		r2 = r * r;

		for (int p = 0; p < spf.numOfSplineParts; p++)
		{
			//first derivative
			spftmp1[spf.numOfSplineParts - 1 - p] = spf.splineWeights[bin - p][p][1] + 2.0 * spf.splineWeights[bin - p][p][2] * r + 3.0 * spf.splineWeights[bin - p][p][3] * r2;
			//second derivative
			spftmp2[spf.numOfSplineParts - 1 - p] = 2.0 * spf.splineWeights[bin - p][p][2] + 6.0 * spf.splineWeights[bin - p][p][3] * r;
		}

		evecr[0] = 1.0; //TODO: is this the unit vector? do not calculate again
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
					//INFO: no potential for contact interaction if gamma is set
					if (gamma == 0.0)
					{
						//square well
						if (rni < potentialRange)
						{
							potentialIntern += potentialStrength;
						}
					}
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

vector<double> InhContactBosons::GetCenterOfMass(vector<vector<double> >& R)
{
	vector<double> com = { 0.0, 0.0, 0.0 };
	return com;
}

double InhContactBosons::GetExternalPotential(vector<double>& r)
{
	double x;
	double k;
	double value = 0.0;
	double potentialK = params[2]; //given in units of kf
	double potentialStrength = params[3]; //given in units of recoil energy E_r
	double kf = M_PI; //INFO: kf = \pi * n = \pi * N / LBOX; change this if density n is not equal to one.
	x = GetCoordinateNIC(r[0]) + halfLength;
	if (params.size() > 4)
	{
		potentialK = params[7];
		potentialStrength = params[8];
	}
	if (potentialK > 0.0)
	{
		k = potentialK * kf;
		if (potentialStrength > 0.0)
		{
			value = sin(k * x);
			value *= value; //sin^2
			value *= potentialStrength; //corresponds to V_0 / E_r
			value *= k * k; //INFO: potential is given by k^2 V_0 / E_r sin^2(k x) to get a factor of 1 in front of kinetic energy
		}
	}
	//INFO: add perturbation
	if (params.size() > 8)
	{
		//additive noise
		//value += params[9] * randomNormal();

		//multiply external potential by gauss in time
		//value *= 1.0 + exp(-pow((this->time - params[9])/(params[9] * 0.4), 2.0)) * 0.001;

		//add external potential of different k and gaussian in time
		double pulseK = params[10];
		double pulseStrength = params[11];
		//for (int i = -5; i <= 4; i++)
		//for (double i = 1; i <= LBOX_2; i += 1)
		{
			//pulseK = params[10] + i * 0.05;
			//pulseK = 1.0 / LBOX * kf; //pulse lowest possible k
			pulseK = 0.5 * kf; //pulse on center of BZ
			double pulse = sin(pulseK * x);
			pulse *= pulse; //sin^2
			pulse -= 0.5; //symmetric around zero
			pulse *= pulseStrength;
			pulse *= k * k;
			pulse *= exp(-pow((this->time - params[9]) / (params[9] * 0.4), 2.0));
			value += pulse;
		}
	}
	if (params.size() > 4 && params[4] > 0.0 && this->time >= params[4])
	{
		value = 0.0;
	}
	else if (params.size() > 8)
	{
		//value *= sin(2.0 * M_PI * this->time * params[9]);
		//value *= (sin(2.0 * M_PI * this->time * params[9]) + sin(2.0 * M_PI * this->time * params[9] * 0.5)) * 0.5;
		//gauss
		//value *= exp(-pow((this->time - params[9])/(params[9] * 0.4), 2.0));
	}
	//value = 0.0;
	//value = abs(x);
	//value = 1.0;
	return value;
}

void InhContactBosons::CalculateExpectationValues(vector<double>& O, vector<vector<vector<double> > >& spf_sD, vector<vector<double> >& spf_sD2, vector<vector<vector<double> > >& pc_sD, vector<vector<double> >& pc_sD2, vector<double>& otherO, vector<double>& uR, vector<double>& uI, double phiR, double phiI)
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
			for (int a = 0; a < DIM; a++)
			{
				//INFO: use the last three splineSumsD with the first three parameters
				tmp = (spf.bcFactorsEnd[i][0] * spf_sD[spf.numberOfSplines - 3][n][a] + spf.bcFactorsEnd[i][1] * spf_sD[spf.numberOfSplines - 2][n][a] + spf.bcFactorsEnd[i][2] * spf_sD[spf.numberOfSplines - 1][n][a]);
				vecKineticSumR1[a] += uR[i] * tmp;
				vecKineticSumI1[a] += uI[i] * tmp;
			}
			//INFO: use the last three splineSumsD2 with the first three parameters
			tmp = (spf.bcFactorsEnd[i][0] * spf_sD2[spf.numberOfSplines - 3][n] + spf.bcFactorsEnd[i][1] * spf_sD2[spf.numberOfSplines - 2][n] + spf.bcFactorsEnd[i][2] * spf_sD2[spf.numberOfSplines - 1][n]);
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
		for (int a = 0; a < DIM; a++)
		{
			tmp = -2.0 * gamma * pc.nodeSpacing * pc_sD[0][n][a];
			vecKineticSumR1[a] += tmp;
			//vecKineticSumI1[a] += tmp;
		}
		tmp = -2.0 * gamma * pc.nodeSpacing * pc_sD2[0][n];
		kineticSumR2 += tmp;
		//kineticSumI2 += tmp;
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

	//cout << "potentialExtern=" << otherO[2] << endl;
	//cout << "potentialIntern=" << otherO[0] << endl;
	//cout << "kineticSumR1=" << kineticSumR1 << endl;
	//cout << "kineticSumI1=" << kineticSumI1 << endl;
	//cout << "kineticSumR2=" << kineticSumR2 << endl;
	//cout << "kineticSumI2=" << kineticSumI2 << endl;
	//cout << "kineticSumR1I1=" << kineticSumR1I1 << endl;
	//cout << "kineticR=" << kineticR << endl;
	//cout << "kineticI=" << kineticI << endl;

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

void InhContactBosons::CalculateExpectationValues(vector<vector<double> >& R, vector<double>& uR, vector<double>& uI, double phiR, double phiI)
{
	RefreshLocalOperators();
	CalculateOtherLocalOperators(R);
	CalculateExpectationValues(this->localOperators, this->spf.splineSumsD, this->spf.splineSumsD2, this->pc.splineSumsD, this->pc.splineSumsD2, this->otherLocalOperators, uR, uI, phiR, phiI);
}

void InhContactBosons::CalculateExpectationValues(ICorrelatedSamplingData* sample, vector<double>& uR, vector<double>& uI, double phiR, double phiI)
{
	cout << "not implemented" << endl;
}

void InhContactBosons::CalculateAdditionalSystemProperties(vector<vector<double> >& R, vector<double>& uR, vector<double>& uI, double phiR, double phiI)
{
	additionalObservables.ClearValues();
	double r;
	double rni;
	double rn;
	double ri;
	vector<double> vecrni(DIM);
	double weight = 1.0;
	double weight2 = 1.0;

	//density
	weight = 1.0 / 2.0;
	for (int i = 0; i < N; i++)
	{
		r = GetCoordinateNIC(R[i][0]) + halfLength; //to get coordinate in interval [0, LBOX]
		if (r < density.grid.max)
		{
			density.AddToHistogram(0, r, weight); //INFO: only for 1D - otherwise there should be a multi-dim histogram variable r
		}
	}

	//pairDensity, pairDistribution
	weight = 1.0 / ((double) (N - 1));
	weight2 = 1.0;
	for (int i = 0; i < N; i++)
	{
		ri = GetCoordinateNIC(R[i][0]) + halfLength; //to get coordinate in interval [0, LBOX]
		for (int j = 0; j < i; j++)
		{
			rn = GetCoordinateNIC(R[j][0]) + halfLength; //to get coordinate in interval [0, LBOX]
			if (ri < pairDensity.grids[0].max && rn < pairDensity.grids[1].max)
			{
				pairDensity.AddToHistogram(0, { ri, rn }, weight2);
			}
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

void InhContactBosons::CalculateWavefunction(vector<double>& O, vector<double>& uR, vector<double>& uI, double phiR, double phiI)
{
	double sum = 0;

	for (int i = 0; i < N_PARAM; i++)
	{
		sum += uR[i] * O[i];
	}
	sum += -2.0 * gamma * pc.nodeSpacing * pc.splineSums[0];

	//cout << scientific << setprecision(16) << "sum:" << sum << endl;
	exponent = sum;
	wf = exp(exponent + phiR);
}

void InhContactBosons::CalculateWavefunction(vector<vector<double> >& R, vector<double>& uR, vector<double>& uI, double phiR, double phiI)
{
	CalculateLocalOperators(R);
	CalculateWavefunction(this->localOperators, uR, uI, phiR, phiI);
}

void InhContactBosons::CalculateWavefunction(ICorrelatedSamplingData* sample, vector<double>& uR, vector<double>& uI, double phiR, double phiI)
{
	CSDataBulkSplines* s = dynamic_cast<CSDataBulkSplines*>(sample);
	CalculateWavefunction(s->localOperators, uR, uI, phiR, phiI);
}

void InhContactBosons::CalculateWFChange(vector<vector<double> >& R, vector<double>& uR, vector<double>& uI, double phiR, double phiI, int changedParticleIndex, vector<double>& oldPosition)
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
	r = GetCoordinateNIC(oldPosition[0]) + halfLength; //to get coordinate in interval [0, LBOX]
	bin = GetBinIndex(spf.nodes, r);
	r2 = r * r;
	r3 = r2 * r;
	for (int p = 0; p < spf.numOfSplineParts; p++)
	{
		spf.sumOldPerBin[bin - p] += spf.splineWeights[bin - p][p][0] + spf.splineWeights[bin - p][p][1] * r + spf.splineWeights[bin - p][p][2] * r2 + spf.splineWeights[bin - p][p][3] * r3;
	}
	r = GetCoordinateNIC(R[changedParticleIndex][0]) + halfLength; //to get coordinate in interval [0, LBOX]
	bin = GetBinIndex(spf.nodes, r);
	r2 = r * r;
	r3 = r2 * r;
	for (int p = 0; p < spf.numOfSplineParts; p++)
	{
		spf.sumNewPerBin[bin - p] += spf.splineWeights[bin - p][p][0] + spf.splineWeights[bin - p][p][1] * r + spf.splineWeights[bin - p][p][2] * r2 + spf.splineWeights[bin - p][p][3] * r3;
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
		sum += uR[i] * (spf.bcFactorsEnd[i][0] * spf.splineSumsNew[spf.numberOfSplines - 3] + spf.bcFactorsEnd[i][1] * spf.splineSumsNew[spf.numberOfSplines - 2] + spf.bcFactorsEnd[i][2] * spf.splineSumsNew[spf.numberOfSplines - 1]);
	}
	spIdx = 3;
	for (int i = spf.np1; i < spf.np2; i++)
	{
		sum += uR[i] * spf.splineSumsNew[spIdx];
		spIdx++;
	}

	//INFO: PairCorrelation
	offset = spf.np3;
	for (int i = 0; i < pc.np1; i++)
	{
		sum += uR[offset + i] * (pc.bcFactorsStart[i][0] * pc.splineSumsNew[0] + pc.bcFactorsStart[i][1] * pc.splineSumsNew[1] + pc.bcFactorsStart[i][2] * pc.splineSumsNew[2]);
	}
	sum += -2.0 * gamma * pc.nodeSpacing * pc.splineSumsNew[0];
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

double InhContactBosons::CalculateWFQuotient(vector<vector<double> >& R, vector<double>& uR, vector<double>& uI, double phiR, double phiI, int changedParticleIndex, vector<double>& oldPosition)
{
	CalculateWFChange(R, uR, uI, phiR, phiI, changedParticleIndex, oldPosition);
	double wfQuotient = exp(2.0 * (exponentNew - exponent));

	//cout << "qu=" << wfQuotient << " (" << exponentNew << " - " << exponent << ")" << endl;

	return wfQuotient;
}

void InhContactBosons::AcceptMove()
{
	wf = wfNew;
	exponent = exponentNew;
	spf.splineSums = spf.splineSumsNew;
	pc.splineSums = pc.splineSumsNew;
}

void InhContactBosons::InitCorrelatedSamplingData(vector<ICorrelatedSamplingData*>& data)
{
}

void InhContactBosons::FillCorrelatedSamplingData(ICorrelatedSamplingData* data)
{
}

}
