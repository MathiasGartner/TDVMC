#include "Bosons1DMixture.h"

#include "../Potentials/Gauss.h"
#include "../Potentials/None.h"
#include "../Potentials/SquareWell.h"
#include "../SplineFactory.h"

namespace PhysicalSystems
{

Bosons1DMixture::Bosons1DMixture(vector<double>& params, string configDirectory) :
		IPhysicalSystem(params, configDirectory)
{
	this->USE_NORMALIZATION_AND_PHASE = true;
	this->USE_NIC = true;
	this->USE_MOVE_COM_TO_ZERO = false;

	halfLength = 0;
	maxDistance = 0;

	numOfOtherLocalOperators = 0;
	numOfPairDistributionValues = 0;
	numOfkValues = 0;

	exponentNew = 0;
}

Bosons1DMixture::~Bosons1DMixture()
{
	for (auto& ppp : this->particlePairProperties)
	{
		delete ppp.potential;
	}
}

void Bosons1DMixture::ExtendNodes(vector<double>& n)
{
	int periodicNodeCount = 3;
	for (int i = 0; i < periodicNodeCount; i++)
	{
		//NEWTODO: check if this is okay (compare to Bosons1D::SetNodes)
		n.insert(n.begin(), -n[2 * i + 1]);
		n.push_back(2.0 * n[n.size() - 1 - i] - n[n.size() - 2 * (i + 1)]);
	}
}

void Bosons1DMixture::SetNodes(vector<double> n)
{
	this->globalNodes = n;
	ExtendNodes(this->globalNodes);
}

void Bosons1DMixture::SetPairDistributionBinCount(double n)
{
	this->numOfPairDistributionValues = n;
}

void Bosons1DMixture::SetParticleType(vector<int> p)
{
	vector<ParticleType> pt;
	int numOfTypeSpecifiedParticles = 0;
	for (auto i : p)
	{
		pt.push_back(static_cast<ParticleType>(i));
		numOfTypeSpecifiedParticles++;
	}
	//INFO: set type of remaining particles to be the last specified type
	for (int i = numOfTypeSpecifiedParticles; i < N; i++)
	{
		pt.push_back(pt.back());
	}
	SetParticleType(pt);
}

void Bosons1DMixture::SetParticleType(vector<ParticleType> p)
{
	this->originalParticleTypes = p;
	InitVector(particleTypes, N, -1);
	InitVector(correlationTypes, N, N, -1);

	InitVector(particleTypeIndexMapping, ParticleType::ParticleType_COUNT, -1);
	int particleIndex = 0;
	for (int i = 0; i < N; i++)
	{
		ParticleType pt = this->originalParticleTypes[i];
		if (particleTypeIndexMapping[pt] == -1)
		{
			particleTypeIndexMapping[pt] = particleIndex;
			particleIndex++;
		}
		this->particleTypes[i] = particleTypeIndexMapping[pt];
	}
	int particleTypeCount = particleIndex;

	InitVector(correlationIndexMapping, ParticleType::ParticleType_COUNT, ParticleType::ParticleType_COUNT, -1);
	int correlationIndex = 0;
	for (int i = 0; i < N; i++)
	{
		int pt1 = this->originalParticleTypes[i];
		for (int j = 0; j < i; j++)
		{
			int pt2 = this->originalParticleTypes[j];
			int cim = correlationIndexMapping[pt1][pt2];
			if (cim == -1)
			{
				cim = correlationIndex;
				correlationIndexMapping[pt1][pt2] = cim;
				correlationIndexMapping[pt2][pt1] = cim;
				PairCorrelation pc;
				pc.FirstParticleType = particleTypeIndexMapping[pt1];
				pc.SecondParticleType = particleTypeIndexMapping[pt2];
				this->pcs.push_back(pc);
				this->correlationNames.push_back(this->ParticleShortNames[pt1] + "-" + this->ParticleShortNames[pt2]);
				correlationIndex++;
			}
			this->correlationTypes[i][j] = cim;
			this->correlationTypes[j][i] = cim;
			this->pcs[cim].numOfParticlePairs++;
		}
	}
	int correlationsCount = correlationIndex;

	this->particleProperties.resize(particleTypeCount);
	this->particlePairProperties.resize(correlationsCount);
	this->oneParticleData.resize(N);
}

void Bosons1DMixture::InitSystem()
{
	int index;

	double tmphbar = 1.0;
	index = this->particleTypeIndexMapping[ParticleType::Boson];
	if (index > -1)
	{
		ParticleProperties& pp = this->particleProperties[index];
		pp.mass = 0.5;
		pp.hbarOver2m = tmphbar * tmphbar / (2.0 * pp.mass);
		pp.kineticEnergyFactor = 1.0;
	}
	index = this->particleTypeIndexMapping[ParticleType::Impurity];
	if (index > -1)
	{
		ParticleProperties& pp = this->particleProperties[index];
		pp.mass = 0.5; //*2.0;
		pp.hbarOver2m = tmphbar * tmphbar / (2.0 * pp.mass);
		pp.kineticEnergyFactor = 87.0 / 41.0; //INFO: mass ratio of m_B and m_I according to Catani 10.1103/PhysRevA.85.023623
	}
	index = this->particleTypeIndexMapping[ParticleType::Dummy];
	if (index > -1)
	{
		ParticleProperties& pp = this->particleProperties[index];
		pp.mass = 0.5;
		pp.hbarOver2m = tmphbar * tmphbar / (2.0 * pp.mass);
		pp.kineticEnergyFactor = 1.0;
	}

	//INFO: Pair potentials
	double defaultStrength = 10.0;
	double defaultRange = 0.1;
	double defaultGamma = 2.0 * 2.0 / 144.0;
	double eta = 10.0 * 0.0 + 1.0;
	defaultGamma = params[1];
	eta = params[0];
	index = this->correlationIndexMapping[ParticleType::Boson][ParticleType::Boson];
	if (index > -1)
	{
		ParticlePairProperties& ppp = this->particlePairProperties[index];
		//ppp.potential = new Potentials::Gauss(defaultStrength, defaultRange);
		ppp.potential = new Potentials::None();
		PairCorrelation& pc = this->pcs[index];
		pc.useContactInteractionBoundaryCondition = true; //params[0] == 0;
		pc.useContactInteractionBoundaryConditionExp = false; //params[0] == 1;
		//pc.contactLengthscale = 0.01;
		pc.gamma = defaultGamma;
	}
	index = this->correlationIndexMapping[ParticleType::Boson][ParticleType::Impurity];
	if (index > -1)
	{
		ParticlePairProperties& ppp = this->particlePairProperties[index];
		//ppp.potential = new Potentials::Gauss(defaultStrength*2.0, defaultRange);
		ppp.potential = new Potentials::None();
		PairCorrelation& pc = this->pcs[index];
		pc.useContactInteractionBoundaryCondition = true;
		pc.gamma = defaultGamma * eta;
	}
	index = this->correlationIndexMapping[ParticleType::Boson][ParticleType::Dummy];
	if (index > -1)
	{
		ParticlePairProperties& ppp = this->particlePairProperties[index];
		ppp.potential = new Potentials::Gauss(defaultStrength, defaultRange);
	}
	index = this->correlationIndexMapping[ParticleType::Impurity][ParticleType::Impurity];
	if (index > -1)
	{
		ParticlePairProperties& ppp = this->particlePairProperties[index];
		//ppp.potential = new Potentials::Gauss(defaultStrength, defaultRange);
		ppp.potential = new Potentials::None();
		PairCorrelation& pc = this->pcs[index];
		pc.useContactInteractionBoundaryCondition = true;
		pc.gamma = defaultGamma;
	}
	index = this->correlationIndexMapping[ParticleType::Impurity][ParticleType::Dummy];
	if (index > -1)
	{
		ParticlePairProperties& ppp = this->particlePairProperties[index];
		ppp.potential = new Potentials::Gauss(defaultStrength, defaultRange);
	}
	index = this->correlationIndexMapping[ParticleType::Dummy][ParticleType::Dummy];
	if (index > -1)
	{
		ParticlePairProperties& ppp = this->particlePairProperties[index];
		ppp.potential = new Potentials::Gauss(defaultStrength, defaultRange);
	}

	numOfOtherLocalOperators = 4; //INFO: potentialIntern, potentialInternComplex, potentialExtern, potentialExternComplex
	otherLocalOperators.resize(numOfOtherLocalOperators);
	numOfkValues = 150;
	numOfOtherExpectationValues = 9;
	numOfAdditionalSystemProperties = numOfOtherExpectationValues;

	halfLength = LBOX / 2.0;

	int numOfSplinedFuncs = this->pcs.size();
	if (globalNodes.empty())
	{
		//INFO: BC
		//this is for SetBoundaryConditions3_1D_OR_2 and SetBoundaryConditions3_1D_CO_2
		for (int i = -3; i < (N_PARAM / numOfSplinedFuncs) + 3; i++)
		{
			globalNodes.push_back((i * halfLength) / ((double) (N_PARAM / numOfSplinedFuncs) - 1));
		}
	}
	maxDistance = globalNodes[globalNodes.size() - 4];

	int requiredParams = 0;
	for (auto& pc : this->pcs)
	{
		pc.numOfParticlePairs_Inv = 1.0 / pc.numOfParticlePairs;

		pc.nodes = this->globalNodes;
		pc.nodeSpacing = pc.nodes[1] - pc.nodes[0];
		//numberOfSplines = N_PARAM + 1; //BC3.1D.CO.1
		//numberOfSplines = N_PARAM + 2; //BC3.1D.CO.2
		//numberOfSplines = N_PARAM + 3; //BC3.1D.CO.2 + BC3.1D.OR.1
		pc.numberOfSplines = pc.nodes.size() - 3 - 1; //3rd order splines -(d+1) (d...spline dimension)
		pc.splineWeights = SplineFactory::GetWeights(pc.nodes);
		pc.Init();

		//cut off
		SplineFactory::SetBoundaryConditions3_1D_OR_2(pc.nodes, pc.bcFactorsStart, true);
		SplineFactory::SetBoundaryConditions3_1D_CO_2(pc.nodes, pc.bcFactorsEnd, maxDistance, true);
		pc.numberOfSpecialParametersStart = pc.bcFactorsStart.size();
		pc.numberOfSpecialParametersEnd = pc.bcFactorsEnd.size();
		pc.numberOfStandardParameters = (N_PARAM / numOfSplinedFuncs) - pc.numberOfSpecialParametersStart - pc.numberOfSpecialParametersEnd;
		pc.numberOfTotalParameters = pc.numberOfSpecialParametersStart + pc.numberOfSpecialParametersEnd + pc.numberOfStandardParameters;
		pc.np1 = pc.numberOfSpecialParametersStart;
		pc.np2 = pc.np1 + pc.numberOfStandardParameters;
		pc.np3 = pc.np2 + pc.numberOfSpecialParametersEnd;

		requiredParams += pc.numberOfSplines - (3 - pc.numberOfSpecialParametersStart) - (3 - pc.numberOfSpecialParametersEnd);
	}

	if (N_PARAM != requiredParams)
	{
		cout << "!!! WRONG NUMBER OF PARAMETERS !!!" << endl;
		cout << "!!! need " << requiredParams << " parameters, got " << N_PARAM << " !!!" << endl;
		exit(0);
	}

	wf = 0.0;
	exponent = 0.0;

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

	//NEWTODO: create observables for pair distribution of bosons and impurities seperately
	pairDistribution.name = "pairDistribution";
	pairDistribution.grid.name = "r_ij";
	//pairDistribution.InitGrid(0.0, halfLength, 0.05);
	pairDistribution.InitGrid(0.0, globalNodes[globalNodes.size() - 4], globalNodes[globalNodes.size() - 4] / numOfPairDistributionValues);
	pairDistribution.InitScaling();
	pairDistribution.InitObservables( { "g_2(r_ij)" });

	structureFactor.name = "structureFactor";
	structureFactor.grid.name = "k";
	structureFactor.InitGrid(kNorms);
	structureFactor.InitObservables( { "S(k)" });

	density.name = "density";
	density.InitGrid(0.0, LBOX, 0.1);
	vector<string> particleNames;
	for (int i = 0; i < this->ParticleType::ParticleType_COUNT; i++)
	{
		int idx = this->particleTypeIndexMapping[i];
		if (idx >= 0)
		{
			particleNames.push_back(this->ParticleNames[idx]);
		}
	}
	density.InitObservables(particleNames);

	pairDistributionsByParticleTypes.name = "pairDistributionsByParticleTypes";
	pairDistributionsByParticleTypes.InitGrid(0.0, globalNodes[globalNodes.size() - 4], globalNodes[globalNodes.size() - 4] / numOfPairDistributionValues);
	pairDistributionsByParticleTypes.InitObservables(this->correlationNames);

	additionalObservables.Add(&pairDistribution);
	additionalObservables.Add(&structureFactor);
	additionalObservables.Add(&density);
	additionalObservables.Add(&pairDistributionsByParticleTypes);
}

void Bosons1DMixture::RefreshLocalOperators()
{
	int offset = 0;
	for (auto& pc : this->pcs)
	{
		//INFO: BC
		for (int i = 0; i < pc.np1; i++)
		{
			this->localOperators[i + offset] = (pc.bcFactorsStart[i][0] * pc.splineSums[0] + pc.bcFactorsStart[i][1] * pc.splineSums[1] + pc.bcFactorsStart[i][2] * pc.splineSums[2]);
		}
		int spIdx = 3;
		for (int i = pc.np1; i < pc.np2; i++)
		{
			this->localOperators[i + offset] = pc.splineSums[spIdx];
			spIdx++;
		}
		int idx = 0;
		for (int i = pc.np2; i < pc.np3; i++)
		{
			this->localOperators[i + offset] = (pc.bcFactorsEnd[idx][0] * pc.splineSums[pc.numberOfSplines - 3] + pc.bcFactorsEnd[idx][1] * pc.splineSums[pc.numberOfSplines - 2] + pc.bcFactorsEnd[idx][2] * pc.splineSums[pc.numberOfSplines - 1]);
			idx++;
		}
		offset += pc.numberOfTotalParameters;
	}
}

void Bosons1DMixture::CalculateLocalOperators(vector<vector<double> >& R)
{
	int ct;

	int bin;
	double rni;
	double rni2;
	double rni3;
	vector<double> vecrni(DIM);

	for (auto& pc : this->pcs)
	{
		ClearVector(pc.splineSums);
	}

	for (int n = 0; n < N; n++)
	{
		for (int i = 0; i < n; i++)
		{
			ct = this->correlationTypes[n][i];
			PairCorrelation& pc = this->pcs[ct];
			rni = VectorDisplacementNIC_DIM(R[n], R[i], vecrni);
			if (rni <= maxDistance)
			{
				bin = pc.BinIndex(rni);
				rni2 = rni * rni;
				rni3 = rni2 * rni;
				for (int p = 0; p < 4; p++)
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

void Bosons1DMixture::CalculateOtherLocalOperators(vector<vector<double> >& R)
{
	int ct; //mapped CorrelationType

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
	/*
	 double a = params[0];
	 double b = params[1];
	 if (params.size() > 2 && this->time >= params[2])
	 {
	 a = params[3];
	 b = params[4];
	 }
	 */

	for (auto& pc : this->pcs)
	{
		ClearVector(pc.splineSumsD);
		ClearVector(pc.splineSumsD2);
		ClearVector(pc.contactSumD);
		ClearVector(pc.contactSumD2);
	}
	for (auto& opd : this->oneParticleData)
	{
		ClearVector(opd.vecKineticSumR1);
		ClearVector(opd.vecKineticSumI1);
		opd.kineticSumR1 = 0;
		opd.kineticSumI1 = 0;
		opd.kineticSumR1I1 = 0;
		opd.kineticSumR2 = 0;
		opd.kineticSumI2 = 0;
	}
	ClearVector(otherLocalOperators);

	for (int n = 0; n < N; n++)
	{
		//external potential energy
		potentialExtern += GetExternalPotential(R[n]);

		for (int i = 0; i < N; i++)
		{
			ct = this->correlationTypes[n][i];
			PairCorrelation& pc = this->pcs[ct];
			ParticlePairProperties& ppp = this->particlePairProperties[ct];

			rni = VectorDisplacementNIC_DIM(R[n], R[i], vecrni);

			if (rni <= maxDistance)
			{
				//internal potential energy
				if (i < n)
				{
					potentialIntern += ppp.potential->GetPotential(rni);
				}

				//kinetic energy
				//TODO: maybe it is faster to calculate all localOperators[i][n] before and then loop over this precalculated values
				//		otherwise values are calculated multiple times
				if (i != n) //TODO: also include in (i < n) branch -> then only loop over i < n in "for (int i = 0; i < N; i++)"
				{
					bin = pc.BinIndex(rni);
					rni2 = rni * rni;

					for (int a = 0; a < DIM; a++)
					{
						evecrni[a] = vecrni[a] / rni;
					}
					bool outsideL2 = false;
					double firstDerivativeFactor = outsideL2 ? -1.0 : 1.0;

					for (int a = 0; a < DIM; a++)
					{
						pc.contactSumD[n][a] += evecrni[a] * pc.gamma * exp(-rni / pc.contactLengthscale); //in 2D or 3D one should use vecrni[a] instead of rni (which could also be negative) - check this!!
					}
					pc.contactSumD2[n] += -pc.gamma / pc.contactLengthscale * exp(-rni / pc.contactLengthscale);

					for (int p = 0; p < s; p++)
					{
						//first derivative
						tmp1[3 - p] = pc.splineWeights[bin - p][p][1] + 2.0 * pc.splineWeights[bin - p][p][2] * rni + 3.0 * pc.splineWeights[bin - p][p][3] * rni2;
						//second derivative
						tmp2[3 - p] = 2.0 * pc.splineWeights[bin - p][p][2] + 6.0 * pc.splineWeights[bin - p][p][3] * rni;
					}
					for (int a = 0; a < DIM; a++)
					{
						for (int b = 0; b < s; b++)
						{
							pc.splineSumsD[bin - b][n][a] += tmp1[3 - b] * evecrni[a] * firstDerivativeFactor;
						}
					}
					double secondDerivativeFactor = DIM - 1.0; //INFO: at least correct for 2D and 3D (and 1D?)
					for (int b = 0; b < s; b++)
					{
						pc.splineSumsD2[bin - b][n] += tmp2[3 - b] + secondDerivativeFactor / rni * tmp1[3 - b];
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

vector<double> Bosons1DMixture::GetCenterOfMass(vector<vector<double> >& R)
{
	vector<double> com = { 0.0, 0.0, 0.0 };
	//TODO: implement
	return com;
}

double Bosons1DMixture::GetExternalPotential(vector<double>& r)
{
	return 0;
}

void Bosons1DMixture::CalculateExpectationValues(vector<double>& O, vector<double>& otherO, vector<double>& gr, vector<double>& uR, vector<double>& uI, double phiR, double phiI)
{
	int pt; //mapped ParticleType

	double tmp = 0;
	double kineticR = 0;
	double kineticI = 0;

	double kineticSumR1SumOPD = 0;
	double kineticSumR2SumOPD = 0;

	localEnergyR = 0;
	localEnergyI = 0;

	ClearVector(localOperatorsMatrix);
	ClearVector(localOperatorlocalEnergyR);
	ClearVector(localOperatorlocalEnergyI);
	ClearVector(otherExpectationValues);

	for (int n = 0; n < N; n++)
	{
		pt = this->particleTypes[n];
		ParticleProperties& pp = this->particleProperties[pt];
		OneParticleData& opd = this->oneParticleData[n];
		for (int a = 0; a < DIM; a++)
		{
			opd.vecKineticSumR1[a] = 0.0;
			opd.vecKineticSumI1[a] = 0.0;
		}

		int offset = 0;
		for (auto& pc : this->pcs)
		{
			bool pcIsRelevant = pc.FirstParticleType == pt || pc.SecondParticleType == pt;
			if (pcIsRelevant)
			{
				//INFO: BC
				if (pc.useContactInteractionBoundaryCondition)
				{
					for (int a = 0; a < DIM; a++)
					{
						tmp = -2.0 * pc.gamma * pc.nodeSpacing * pc.splineSumsD[0][n][a];
						opd.vecKineticSumR1[a] += tmp;
						//opd.vecKineticSumI1[a] += tmp; //INFO: imaginary part of boundary condition is zero
					}
					tmp = -2.0 * pc.gamma * pc.nodeSpacing * pc.splineSumsD2[0][n];
					opd.kineticSumR2 += tmp;
					//kineticSumI2 += tmp; //INFO: imaginary part of boundary condition is zero
				}
				else if (pc.useContactInteractionBoundaryConditionExp)
				{
					for (int a = 0; a < DIM; a++)
					{
						tmp = 2.0 * pc.gamma * exp(-LBOX_2 / pc.contactLengthscale) * pc.nodeSpacing * pc.splineSumsD.back()[n][a];
						opd.vecKineticSumR1[a] += tmp;
						tmp = pc.contactSumD[n][a];
						opd.vecKineticSumR1[a] += tmp;
						//opd.vecKineticSumI1[a] += tmp; //INFO: imaginary part of boundary condition is zero
					}
					tmp = 2.0 * pc.gamma * exp(-LBOX_2 / pc.contactLengthscale) * pc.nodeSpacing * pc.splineSumsD2.back()[n];
					opd.kineticSumR2 += tmp;
					tmp = pc.contactSumD2[n];
					opd.kineticSumR2 += tmp;
					//kineticSumI2 += tmp; //INFO: imaginary part of boundary condition is zero
				}
				for (int i = 0; i < pc.np1; i++)
				{
					for (int a = 0; a < DIM; a++)
					{
						tmp = (pc.bcFactorsStart[i][0] * pc.splineSumsD[0][n][a] + pc.bcFactorsStart[i][1] * pc.splineSumsD[1][n][a] + pc.bcFactorsStart[i][2] * pc.splineSumsD[2][n][a]);
						opd.vecKineticSumR1[a] += uR[i + offset] * tmp;
						opd.vecKineticSumI1[a] += uI[i + offset] * tmp;
					}
					tmp = (pc.bcFactorsStart[i][0] * pc.splineSumsD2[0][n] + pc.bcFactorsStart[i][1] * pc.splineSumsD2[1][n] + pc.bcFactorsStart[i][2] * pc.splineSumsD2[2][n]);
					opd.kineticSumR2 += uR[i + offset] * tmp;
					opd.kineticSumI2 += uI[i + offset] * tmp;
				}
				int spIdx = 3;
				for (int i = pc.np1; i < pc.np2; i++)
				{
					for (int a = 0; a < DIM; a++)
					{
						opd.vecKineticSumR1[a] += uR[i + offset] * pc.splineSumsD[spIdx][n][a];
						opd.vecKineticSumI1[a] += uI[i + offset] * pc.splineSumsD[spIdx][n][a];
					}
					opd.kineticSumR2 += uR[i + offset] * pc.splineSumsD2[spIdx][n];
					opd.kineticSumI2 += uI[i + offset] * pc.splineSumsD2[spIdx][n];
					spIdx++;
				}
				int idx = 0;
				for (int i = pc.np2; i < pc.np3; i++)
				{
					for (int a = 0; a < DIM; a++)
					{
						tmp = (pc.bcFactorsEnd[idx][0] * pc.splineSumsD[pc.numberOfSplines - 3][n][a] + pc.bcFactorsEnd[idx][1] * pc.splineSumsD[pc.numberOfSplines - 2][n][a] + pc.bcFactorsEnd[idx][2] * pc.splineSumsD[pc.numberOfSplines - 1][n][a]);
						opd.vecKineticSumR1[a] += uR[i + offset] * tmp;
						opd.vecKineticSumI1[a] += uI[i + offset] * tmp;
					}
					tmp = (pc.bcFactorsEnd[idx][0] * pc.splineSumsD2[pc.numberOfSplines - 3][n] + pc.bcFactorsEnd[idx][1] * pc.splineSumsD2[pc.numberOfSplines - 2][n] + pc.bcFactorsEnd[idx][2] * pc.splineSumsD2[pc.numberOfSplines - 1][n]);
					opd.kineticSumR2 += uR[i + offset] * tmp;
					opd.kineticSumI2 += uI[i + offset] * tmp;
					idx++;
				}
			}
			offset += pc.numberOfTotalParameters;
		}

		opd.kineticSumR1I1 += 2.0 * VectorDotProduct_DIM(opd.vecKineticSumR1, opd.vecKineticSumI1);
		opd.kineticSumR1 += VectorNorm2_DIM(opd.vecKineticSumR1);
		opd.kineticSumI1 += VectorNorm2_DIM(opd.vecKineticSumI1);

		kineticR += -pp.kineticEnergyFactor * (opd.kineticSumR1 - opd.kineticSumI1 + opd.kineticSumR2);
		kineticI += -pp.kineticEnergyFactor * (opd.kineticSumR1I1 + opd.kineticSumI2);

		kineticSumR1SumOPD += opd.kineticSumR1;
		kineticSumR2SumOPD += opd.kineticSumR2;
	}

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
	otherExpectationValues[4] = kineticSumR1SumOPD; //kineticSumR1;
	otherExpectationValues[5] = 0; //kineticSumI1;
	otherExpectationValues[6] = kineticSumR2SumOPD; //kineticSumR2;
	otherExpectationValues[7] = 0; //kineticSumI2;
	otherExpectationValues[8] = 0; //kineticSumR1I1;
}

void Bosons1DMixture::CalculateExpectationValues(vector<vector<double> >& R, vector<double>& uR, vector<double>& uI, double phiR, double phiI)
{
	RefreshLocalOperators();
	CalculateOtherLocalOperators(R);
	vector<double> dummy_grBins;
	CalculateExpectationValues(this->localOperators, this->otherLocalOperators, dummy_grBins, uR, uI, phiR, phiI);
}

void Bosons1DMixture::CalculateExpectationValues(ICorrelatedSamplingData* sample, vector<double>& uR, vector<double>& uI, double phiR, double phiI)
{
	//TODO: implement
}

void Bosons1DMixture::CalculateAdditionalSystemProperties(vector<vector<double> >& R, vector<double>& uR, vector<double>& uI, double phiR, double phiI)
{
	additionalObservables.ClearValues();
	double rni;
	vector<double> vecrni(DIM);

	//pairDistribution
	double weight = 1.0 / ((double) (N - 1));
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < i; j++)
		{
			rni = VectorDisplacementNIC_DIM(R[i], R[j], vecrni);
			if (rni < pairDistribution.grid.max)
			{
				pairDistribution.AddToHistogram(0, rni, weight);
			}

			int ct = this->correlationTypes[i][j];
			if (rni < pairDistributionsByParticleTypes.grid.max)
			{
				pairDistributionsByParticleTypes.AddToHistogram(ct, rni, this->pcs[ct].numOfParticlePairs_Inv);
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

	//density
	for (int i = 0; i < N; i++)
	{
		int pt = this->particleTypes[i];
		double r = GetCoordinateNIC(R[i][0]) + LBOX_2; //to get results in the range [0, L/2]
		density.AddToHistogram(this->particleTypeIndexMapping[pt], r, 1.0);
	}
}

void Bosons1DMixture::CalculateWavefunction(vector<double>& O, vector<double>& uR, vector<double>& uI, double phiR, double phiI)
{
	double sum = 0;

	for (auto& pc : this->pcs)
	{
		if (pc.useContactInteractionBoundaryCondition)
		{
			sum += -2.0 * pc.gamma * pc.nodeSpacing * pc.splineSums[0];
		}
		else if (pc.useContactInteractionBoundaryConditionExp)
		{
			sum += 2.0 * pc.gamma * exp(-LBOX_2 / pc.contactLengthscale) * pc.nodeSpacing * pc.splineSums.back();
			sum += pc.contactSum;
		}
	}
	for (int i = 0; i < N_PARAM; i++)
	{
		sum += uR[i] * O[i];
	}

	//cout << scientific << setprecision(16) << "sum:" << sum << endl;
	exponent = sum;
	wf = exp(exponent + phiR);
}

void Bosons1DMixture::CalculateWavefunction(vector<vector<double> >& R, vector<double>& uR, vector<double>& uI, double phiR, double phiI)
{
	CalculateLocalOperators(R);
	CalculateWavefunction(this->localOperators, uR, uI, phiR, phiI);
}

void Bosons1DMixture::CalculateWavefunction(ICorrelatedSamplingData* sample, vector<double>& uR, vector<double>& uI, double phiR, double phiI)
{
	//TODO: implement
}

void Bosons1DMixture::CalculateWFChange(vector<vector<double> >& R, vector<double>& uR, vector<double>& uI, double phiR, double phiI, int changedParticleIndex, vector<double>& oldPosition)
{
	int ct;

	double sum = 0;
	int bin;
	double rni;
	double rni2;
	double rni3;
	vector<double> vecrni(DIM);

	for (auto& pc : this->pcs)
	{
		ClearVector(pc.sumOldPerBin);
		ClearVector(pc.sumNewPerBin);
		pc.contactSum = 0.0;
		pc.contactOld = 0.0;
		pc.contactNew = 0.0;
	}

	for (int i = 0; i < N; i++)
	{
		if (i != changedParticleIndex)
		{
			ct = this->correlationTypes[i][changedParticleIndex];
			PairCorrelation& pc = this->pcs[ct];
			rni = VectorDisplacementNIC_DIM(R[i], oldPosition, vecrni);
			if (rni <= maxDistance)
			{
				bin = pc.BinIndex(rni);
				if (bin >= pc.numberOfSplines)
				{
					cout << "!!!!" << endl;
				}
				rni2 = rni * rni;
				rni3 = rni2 * rni;
				for (int p = 0; p < 4; p++)
				{
					pc.sumOldPerBin[bin - p] += pc.splineWeights[bin - p][p][0] + pc.splineWeights[bin - p][p][1] * rni + pc.splineWeights[bin - p][p][2] * rni2 + pc.splineWeights[bin - p][p][3] * rni3;
				}
				pc.contactOld += -pc.gamma * pc.contactLengthscale * exp(-rni / pc.contactLengthscale);
			}
			else
			{
				cout << "!!!" << endl;
			}
			rni = VectorDisplacementNIC_DIM(R[i], R[changedParticleIndex], vecrni);
			if (rni <= maxDistance)
			{
				bin = pc.BinIndex(rni);
				rni2 = rni * rni;
				rni3 = rni2 * rni;
				for (int p = 0; p < 4; p++)
				{
					pc.sumNewPerBin[bin - p] += pc.splineWeights[bin - p][p][0] + pc.splineWeights[bin - p][p][1] * rni + pc.splineWeights[bin - p][p][2] * rni2 + pc.splineWeights[bin - p][p][3] * rni3;
				}
				pc.contactNew += -pc.gamma * pc.contactLengthscale * exp(-rni / pc.contactLengthscale);
			}
			else
			{
				cout << "!!!" << endl;
			}
		}
	}

	for (auto& pc : this->pcs)
	{
		for (int i = 0; i < pc.numberOfSplines; i++)
		{
			pc.splineSumsNew[i] = max(0.0, pc.splineSums[i] - pc.sumOldPerBin[i] + pc.sumNewPerBin[i]);
		}
		pc.contactSumNew = min(0.0, pc.contactSum - pc.contactOld + pc.contactNew); //INFO: contact sum hast to be negative, there is no variational parameter as in the case of spline sums
	}

	int offset = 0;
	for (auto& pc : this->pcs)
	{
		//INFO: BCx
		if (pc.useContactInteractionBoundaryCondition)
		{
			sum += -2.0 * pc.gamma * pc.nodeSpacing * pc.splineSumsNew[0];
		}
		else if (pc.useContactInteractionBoundaryConditionExp)
		{
			sum += 2.0 * pc.gamma * exp(-LBOX_2 / pc.contactLengthscale) * pc.nodeSpacing * pc.splineSumsNew.back();
			sum += pc.contactSumNew;
		}
		for (int i = 0; i < pc.np1; i++)
		{
			sum += uR[i + offset] * (pc.bcFactorsStart[i][0] * pc.splineSumsNew[0] + pc.bcFactorsStart[i][1] * pc.splineSumsNew[1] + pc.bcFactorsStart[i][2] * pc.splineSumsNew[2]);
		}
		int spIdx = 3;
		for (int i = pc.np1; i < pc.np2; i++)
		{
			sum += uR[i + offset] * pc.splineSumsNew[spIdx];
			spIdx++;
		}
		int idx = 0;
		for (int i = pc.np2; i < pc.np3; i++)
		{
			sum += uR[i + offset] * (pc.bcFactorsEnd[idx][0] * pc.splineSumsNew[pc.numberOfSplines - 3] + pc.bcFactorsEnd[idx][1] * pc.splineSumsNew[pc.numberOfSplines - 2] + pc.bcFactorsEnd[idx][2] * pc.splineSumsNew[pc.numberOfSplines - 1]);
			idx++;
		}
		offset += pc.numberOfTotalParameters;
	}

	exponentNew = sum;
	wfNew = exp(exponentNew + phiR);
}

double Bosons1DMixture::CalculateWFQuotient(vector<vector<double> >& R, vector<double>& uR, vector<double>& uI, double phiR, double phiI, int changedParticleIndex, vector<double>& oldPosition)
{
	CalculateWFChange(R, uR, uI, phiR, phiI, changedParticleIndex, oldPosition);
	double wfQuotient = exp(2.0 * (exponentNew - exponent));

	//cout << "qu=" << wfQuotient << " (" << exponentNew << " - " << exponent << ")" << endl;

	return wfQuotient;
}

void Bosons1DMixture::AcceptMove()
{
	wf = wfNew;
	exponent = exponentNew;
	for (auto& pc : this->pcs)
	{
		pc.splineSums = pc.splineSumsNew;
		pc.contactSum = pc.contactSumNew;
	}
}

void Bosons1DMixture::InitCorrelatedSamplingData(vector<ICorrelatedSamplingData*>& data)
{
	//TODO: implement
}

void Bosons1DMixture::FillCorrelatedSamplingData(ICorrelatedSamplingData* data)
{
	//TODO: implement
}

}
