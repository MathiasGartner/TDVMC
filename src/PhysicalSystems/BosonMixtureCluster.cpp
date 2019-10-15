#include "BosonMixtureCluster.h"

#include "../Potentials/HFDB_He_He.h"
#include "../Potentials/KTTY_He_Cs.h"
#include "../Potentials/KTTY_He_Na.h"
#include "../Potentials/LJ_He_He.h"
#include "../SplineFactory.h"

using namespace std;

namespace PhysicalSystems
{

BosonMixtureCluster::BosonMixtureCluster(vector<double>& params, string configDirectory) :
		IPhysicalSystem(params, configDirectory)
{
	this->USE_NORMALIZATION_AND_PHASE = true;
	this->USE_NIC = false;
	this->USE_MOVE_COM_TO_ZERO = true;

	exponentNew = 0;

	globalDensityProfileMaxDistance = 0;
	globalNumOfDensityProfileValues = 0;
	globalDensityProfileBinInterval = 0;
	globalDensityProfileBin = 0;
	globalDensityProfileNodePointSpacing = 0;
}

BosonMixtureCluster::~BosonMixtureCluster()
{
	for (auto& ppp : this->particlePairProperties)
	{
		delete ppp.potential;
	}
}

void BosonMixtureCluster::SetNodes(vector<double> n)
{
	this->globalNodes = n;
}

void BosonMixtureCluster::SetDensityProfileBinCount(double n)
{
	this->globalNumOfDensityProfileValues = n;
}

void BosonMixtureCluster::SetParticleType(vector<int> p)
{
	vector<ParticleType> pt;
	for (auto i : p)
	{
		pt.push_back(static_cast<ParticleType>(i));
	}
	SetParticleType(pt);
}

void BosonMixtureCluster::SetParticleType(vector<ParticleType> p)
{
	this->originalParticleTypes = p;
	InitVector(particleTypes, N, 0);
	InitVector(correlationTypes, N, N, 0);

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
			if (correlationIndexMapping[pt1][pt2] == -1)
			{
				correlationIndexMapping[pt1][pt2] = correlationIndex;
				correlationIndexMapping[pt2][pt1] = correlationIndex;
				correlationIndex++;
			}
			this->correlationTypes[i][j] = correlationIndexMapping[pt1][pt2];
			this->correlationTypes[j][i] = correlationIndexMapping[pt2][pt1];
		}
	}
	int correlationsCount = correlationIndex;

	this->corrFuncData.resize(correlationsCount);
	this->particleProperties.resize(particleTypeCount);
	this->particlePairProperties.resize(correlationsCount);
	this->oneParticleData.resize(N);
}

void BosonMixtureCluster::InitSystem()
{
	int index;

	index = this->particleTypeIndexMapping[ParticleType::He3];
	if (index > -1)
	{
		ParticleProperties& pp = this->particleProperties[index];
		pp.mass = 3.0160293191;
		pp.hbarOver2m = pow(Constants::Split::hbar, 2.0) / (2.0 * pp.mass * Constants::Split::u) / (pow(Constants::Split::A2m, 2.0) * Constants::Split::kb);
	}
	index = this->particleTypeIndexMapping[ParticleType::He4];
	if (index > -1)
	{
		ParticleProperties& pp = this->particleProperties[index];
		pp.mass = 4.00260325415;
		pp.hbarOver2m = pow(Constants::Split::hbar, 2.0) / (2.0 * pp.mass * Constants::Split::u) / (pow(Constants::Split::A2m, 2.0) * Constants::Split::kb);
	}
	index = this->particleTypeIndexMapping[ParticleType::Na];
	if (index > -1)
	{
		ParticleProperties& pp = this->particleProperties[index];
		pp.mass = 22.9897692809;
		pp.hbarOver2m = pow(Constants::Split::hbar, 2.0) / (2.0 * pp.mass * Constants::Split::u) / (pow(Constants::Split::A2m, 2.0) * Constants::Split::kb);
	}
	index = this->particleTypeIndexMapping[ParticleType::Cs];
	if (index > -1)
	{
		ParticleProperties& pp = this->particleProperties[index];
		pp.mass = 132.905451932;
		pp.hbarOver2m = pow(Constants::Split::hbar, 2.0) / (2.0 * pp.mass * Constants::Split::u) / (pow(Constants::Split::A2m, 2.0) * Constants::Split::kb);
	}

	//TODO!!!
	index = this->correlationIndexMapping[He4][He4];
	if (index > -1)
	{
		CorrelationFunctionData& cfd = this->corrFuncData[index];
		cfd.mcMillanFactor = -4.7;
		cfd.numberOfSplines = 26; //as in calculations performed for Split X_4 Cluster (BosonClusterWithLogParam)
		cfd.nodes = this->globalNodes; //TODO: allow different nodes
		cfd.nodes -= 0.7;
		cfd.rijSplit = cfd.nodes[3];
		cfd.rijTail = cfd.nodes[cfd.nodes.size() - 4];
		cfd.splineWeights = SplineFactory::GetWeights(cfd.nodes);
		SplineFactory::SetBoundaryConditions1_MM_1(cfd.nodes, cfd.bcFactors, cfd.rijSplit, cfd.mcMillanFactor);
		SplineFactory::SetBoundaryConditions1_EXP_2(cfd.nodes, cfd.bcFactors, cfd.rijTail);
		cfd.Init();
	}

	//INFO: copy
	int copyIndex = this->correlationIndexMapping[He3][He4];
	index = this->correlationIndexMapping[He3][He3];
	if (index > -1)
	{
		this->corrFuncData[index] = this->corrFuncData[copyIndex];
	}
	index = this->correlationIndexMapping[Na][Na];
	if (index > -1)
	{
		this->corrFuncData[index] = this->corrFuncData[copyIndex];
	}
	index = this->correlationIndexMapping[He3][He4];
	if (index > -1)
	{
		CorrelationFunctionData& cfd = this->corrFuncData[index];
		cfd.mcMillanFactor = -4.7;
		cfd.numberOfSplines = 26; //as in calculations performed for Split X_4 Cluster (BosonClusterWithLogParam)
		cfd.nodes = this->globalNodes; //TODO: allow different nodes
		cfd.nodes -= 0.7;
		cfd.rijSplit = cfd.nodes[3];
		cfd.rijTail = cfd.nodes[cfd.nodes.size() - 4];
		cfd.splineWeights = SplineFactory::GetWeights(cfd.nodes);
		SplineFactory::SetBoundaryConditions1_MM_1(cfd.nodes, cfd.bcFactors, cfd.rijSplit, cfd.mcMillanFactor);
		SplineFactory::SetBoundaryConditions1_EXP_2(cfd.nodes, cfd.bcFactors, cfd.rijTail);
		cfd.Init();
	}
	index = this->correlationIndexMapping[He3][Na];
	if (index > -1)
	{
		this->corrFuncData[index] = this->corrFuncData[copyIndex];
	}
	index = this->correlationIndexMapping[He4][Na];
	if (index > -1)
	{
		CorrelationFunctionData& cfd = this->corrFuncData[index];
		cfd.mcMillanFactor = -4.7;
		cfd.numberOfSplines = 26; //as in calculations performed for Split X_4 Cluster (BosonClusterWithLogParam)
		cfd.nodes = this->globalNodes; //TODO: allow different nodes
		cfd.nodes += 1.3;
		cfd.rijSplit = cfd.nodes[3];
		cfd.rijTail = cfd.nodes[cfd.nodes.size() - 4];
		cfd.splineWeights = SplineFactory::GetWeights(cfd.nodes);
		SplineFactory::SetBoundaryConditions1_MM_1(cfd.nodes, cfd.bcFactors, cfd.rijSplit, cfd.mcMillanFactor);
		SplineFactory::SetBoundaryConditions1_EXP_2(cfd.nodes, cfd.bcFactors, cfd.rijTail);
		cfd.Init();
	}
	index = this->correlationIndexMapping[He4][Cs];
	if (index > -1)
	{
		CorrelationFunctionData& cfd = this->corrFuncData[index];
		cfd.mcMillanFactor = -4.7;
		cfd.numberOfSplines = 26; //as in calculations performed for Split X_4 Cluster (BosonClusterWithLogParam)
		cfd.nodes = this->globalNodes; //TODO: allow different nodes
		cfd.nodes += 2.1;
		cfd.rijSplit = cfd.nodes[3];
		cfd.rijTail = cfd.nodes[cfd.nodes.size() - 4];
		cfd.splineWeights = SplineFactory::GetWeights(cfd.nodes);
		SplineFactory::SetBoundaryConditions1_MM_1(cfd.nodes, cfd.bcFactors, cfd.rijSplit, cfd.mcMillanFactor);
		SplineFactory::SetBoundaryConditions1_EXP_2(cfd.nodes, cfd.bcFactors, cfd.rijTail);
		cfd.Init();
	}
	index = this->correlationIndexMapping[He3][Cs];
	if (index > -1)
	{
		CorrelationFunctionData& cfd = this->corrFuncData[index];
		cfd.mcMillanFactor = -4.7;
		cfd.numberOfSplines = 26; //as in calculations performed for Split X_4 Cluster (BosonClusterWithLogParam)
		cfd.nodes = this->globalNodes; //TODO: allow different nodes
		cfd.nodes += 2.1;
		cfd.rijSplit = cfd.nodes[3];
		cfd.rijTail = cfd.nodes[cfd.nodes.size() - 4];
		cfd.splineWeights = SplineFactory::GetWeights(cfd.nodes);
		SplineFactory::SetBoundaryConditions1_MM_1(cfd.nodes, cfd.bcFactors, cfd.rijSplit, cfd.mcMillanFactor);
		SplineFactory::SetBoundaryConditions1_EXP_2(cfd.nodes, cfd.bcFactors, cfd.rijTail);
		cfd.Init();
	}

	//INFO: Pair potentials
	index = this->correlationIndexMapping[He4][He4];
	if (index > -1)
	{
		ParticlePairProperties& ppp = this->particlePairProperties[index];
		ppp.potential = new Potentials::HFDB_He_He();
		//ppp.potential = new Potentials::LJ_He_He();
	}
	index = this->correlationIndexMapping[He3][He3];
	if (index > -1)
	{
		ParticlePairProperties& ppp = this->particlePairProperties[index];
		ppp.potential = new Potentials::HFDB_He_He();
	}
	index = this->correlationIndexMapping[Na][Na];
	if (index > -1)
	{
		//TODO: what potential to use here?
		ParticlePairProperties& ppp = this->particlePairProperties[index];
		ppp.potential = new Potentials::KTTY_He_Na();
	}
	index = this->correlationIndexMapping[He3][He4];
	if (index > -1)
	{
		ParticlePairProperties& ppp = this->particlePairProperties[index];
		ppp.potential = new Potentials::HFDB_He_He();
	}
	index = this->correlationIndexMapping[He3][Na];
	if (index > -1)
	{
		ParticlePairProperties& ppp = this->particlePairProperties[index];
		ppp.potential = new Potentials::KTTY_He_Na();
	}
	index = this->correlationIndexMapping[He4][Na];
	if (index > -1)
	{
		ParticlePairProperties& ppp = this->particlePairProperties[index];
		ppp.potential = new Potentials::KTTY_He_Na();
	}
	index = this->correlationIndexMapping[He4][Cs];
	if (index > -1)
	{
		ParticlePairProperties& ppp = this->particlePairProperties[index];
		ppp.potential = new Potentials::KTTY_He_Cs();
	}
	index = this->correlationIndexMapping[He3][Cs];
	if (index > -1)
	{
		ParticlePairProperties& ppp = this->particlePairProperties[index];
		ppp.potential = new Potentials::KTTY_He_Cs();
	}

	//numOfDensityProfileValues = 200;
	numOfOtherExpectationValues = 3 + globalNumOfDensityProfileValues;

	wf = 0.0;
	exponent = 0.0;

	localEnergyR = 0;
	localEnergyI = 0;
	localOperators.resize(N_PARAM);
	ClearVector (localOperators);
	localOperatorsMatrix.resize(N_PARAM);
	for (auto &row : localOperatorsMatrix)
	{
		row.resize(N_PARAM);
	}
	ClearVector (localOperatorsMatrix);
	localOperatorlocalEnergyR.resize(N_PARAM);
	ClearVector (localOperatorlocalEnergyR);
	localOperatorlocalEnergyI.resize(N_PARAM);
	ClearVector (localOperatorlocalEnergyI);
	otherExpectationValues.resize(numOfOtherExpectationValues);
	ClearVector (otherExpectationValues);

	wfNew = 0.0;
	exponentNew = 0.0;

	globalDensityProfileMaxDistance = globalNodes.back() * 5.0;
	globalDensityProfileNodePointSpacing = globalDensityProfileMaxDistance / globalNumOfDensityProfileValues;
	globalDensityProfileBins.resize(globalNumOfDensityProfileValues);

	globalDensityProfileBinVolumes.resize(globalNumOfDensityProfileValues);
	ClearVector(globalDensityProfileBinVolumes);
	for (int i = 0; i < globalNumOfDensityProfileValues; i++)
	{
		globalDensityProfileBinVolumes[i] = 4.0 * M_PI * pow(globalDensityProfileNodePointSpacing * (i + 1), 3.0) / 3.0;
	}
	for (int i = globalNumOfDensityProfileValues - 1; i > 0; i--)
	{
		globalDensityProfileBinVolumes[i] = globalDensityProfileBinVolumes[i] - globalDensityProfileBinVolumes[i - 1];
	}

	//cout << "nodePointSpacingShort=" << nodePointSpacingShort << ", nodePointSpacingLarge=" << nodePointSpacingLarge << ", rijSplit=" << rijSplit << ", rijSplineSplit=" << rijSplineSplit << ", rijTail=" << rijTail << endl;

	r2.name = "r2";

	angularDistribution.name = "angularDistribution";
	angularDistribution.InitGrid(0.0, 180.0, 1.0);
	angularDistribution.InitObservables({"1-2-3", "1-3-2", "2-1-3"}); //TODO: do not hardcode

	densityFromCOM.name = "densityFromCOM";
	densityFromCOM.InitGrid(0.0, 80.0, 0.1);
	densityFromCOM.InitScaling();
	densityFromCOM.InitObservables({"1", "2", "3"}); //TODO: do not hardcode

	particleDistances.name = "particleDistances";
	particleDistances.InitGrid(0.0, 80, 0.5);
	particleDistances.InitObservables({"1-2", "1-3", "2-3"}); //TODO: do not hardcode

	additionalObservables.Add(&r2);
	additionalObservables.Add(&angularDistribution);
	additionalObservables.Add(&densityFromCOM);
	additionalObservables.Add(&particleDistances);
}

vector<double> BosonMixtureCluster::GetCenterOfMass(vector<vector<double> >& R)
{
	int pt; //mapped ParticleType
	double massSum = 0.0;
	vector<double> com = { 0, 0, 0 };
	for (int i = 0; i < N; i++)
	{
		pt = this->particleTypes[i];
		ParticleProperties& pp = this->particleProperties[pt];
		massSum += pp.mass;
		for (int a = 0; a < DIM; a++)
		{
			com[a] += pp.mass * R[i][a];
		}
	}
	for (int a = 0; a < DIM; a++)
	{
		com[a] /= massSum;
	}
	return com;
}

void BosonMixtureCluster::CalculateOtherLocalOperators(vector<vector<double> >& R)
{
	//TODO: implement
}

void BosonMixtureCluster::CalculateExpectationValues(vector<vector<double> >& R, vector<double>& uR, vector<double>& uI, double phiR, double phiI)
{
	int ct; //mapped CorrelationType
	int pt; //mapped ParticleType

	double temp;

	int s = 4;
	vector<double> vecrni(DIM);
	vector<double> evecrni(DIM);
	vector<double> tmp1(s);
	vector<double> tmp2(s);
	int bin;
	double rni;
	double rni2;

	double potentialExtern = 0;
	double potentialIntern = 0;

	double kineticR = 0;
	double kineticR_part1 = 0;
	double kineticR_part2 = 0;
	double kineticI = 0;

	for (auto& cfd : this->corrFuncData)
	{
		ClearVector(cfd.mcMillanSumD);
		ClearVector(cfd.mcMillanSumD2);
		ClearVector(cfd.constSumD);
		ClearVector(cfd.constSumD2);
		ClearVector(cfd.logSumD);
		ClearVector(cfd.logSumD2);
		ClearVector(cfd.linearSumD);
		ClearVector(cfd.linearSumD2);
		ClearVector(cfd.splineSumsD);
		ClearVector(cfd.splineSumsD2);
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

	localEnergyR = 0;
	localEnergyI = 0;
	ClearVector (localOperators);
	ClearVector (localOperatorsMatrix);
	ClearVector (localOperatorlocalEnergyR);
	ClearVector (localOperatorlocalEnergyI);
	ClearVector (otherExpectationValues);
	ClearVector(globalDensityProfileBins);
	vector<double> com = GetCenterOfMass(R);

	for (int n = 0; n < N; n++)
	{
		pt = this->particleTypes[n];
		ParticleProperties& pp = this->particleProperties[pt];
		OneParticleData& opd = this->oneParticleData[n];

		//external potential energy
		potentialExtern += 0;

		for (int i = 0; i < N; i++)
		{
			ct = this->correlationTypes[n][i];
			CorrelationFunctionData& cfd = this->corrFuncData[ct];
			ParticlePairProperties& ppp = this->particlePairProperties[ct];

			rni = VectorDisplacement(R[n], R[i], vecrni);
			//if (rni < maxDistance)
			{
				//internal potential energy
				if (i < n)
				{
					potentialIntern += ppp.potential->GetPotential(rni);
				}

				//kinetic energy
				//TODO: maybe it is faster to calculate all localOperators[i][n] before and then loop over this precalculated values
				//		otherwise values are calculated multiple times
				if (i != n) //TODO: also include in (i < n) branch -> then only loop over i < n
				{
					for (int a = 0; a < DIM; a++)
					{
						evecrni[a] = vecrni[a] / rni;
					}
					if (rni < cfd.rijSplit)
					{
						double rniMinusD = pow(rni, cfd.mcMillanFactor - 2.0);
						for (int a = 0; a < DIM; a++)
						{
							cfd.mcMillanSumD[n][a] += cfd.mcMillanFactor * rniMinusD * vecrni[a];
						}
						cfd.mcMillanSumD2[n] += cfd.mcMillanFactor * (cfd.mcMillanFactor + 1.0) * rniMinusD;
					}
					else if (rni >= cfd.rijTail)
					{
						for (int a = 0; a < DIM; a++)
						{
							evecrni[a] = vecrni[a] / rni;
						}
						for (int a = 0; a < DIM; a++)
						{
							cfd.constSumD[n][a] += 0.0;
							cfd.linearSumD[n][a] += evecrni[a];
						}
						cfd.constSumD2[n] += 0.0;
						cfd.linearSumD2[n] += 2.0 / rni;
					}
					else
					{
						auto interval = lower_bound(cfd.nodes.begin(), cfd.nodes.end(), rni);
						bin = interval - cfd.nodes.begin() - 1;
						rni2 = rni * rni;

						for (int p = 0; p < s; p++)
						{
							//first derivative
							tmp1[3 - p] = cfd.splineWeights[bin - p][p][1] + 2.0 * cfd.splineWeights[bin - p][p][2] * rni + 3.0 * cfd.splineWeights[bin - p][p][3] * rni2;
							//second derivative
							tmp2[3 - p] = 2.0 * cfd.splineWeights[bin - p][p][2] + 6.0 * cfd.splineWeights[bin - p][p][3] * rni;
						}
						for (int a = 0; a < DIM; a++)
						{
							for (int b = 0; b < s; b++)
							{
								cfd.splineSumsD[bin - b][n][a] += tmp1[3 - b] * evecrni[a];
							}
						}
						for (int b = 0; b < s; b++)
						{
							double secondDerivativeFactor = DIM - 1.0; //INFO: at least correct for 2D and 3D (and 1D?)
							cfd.splineSumsD2[bin - b][n] += tmp2[3 - b] + secondDerivativeFactor / rni * tmp1[3 - b];
						}
					}
					for (int a = 0; a < DIM; a++)
					{
						cfd.logSumD[n][a] += 1.0 / rni * evecrni[a];
					}
					cfd.logSumD2[n] += pow(rni, -2);
				}
			}
		}

		//density profile
		double r;
		vector<double> diff = R[n] - com;
		r = VectorNorm(diff);
		if (r < globalDensityProfileMaxDistance)
		{
			globalDensityProfileBinInterval = r / globalDensityProfileNodePointSpacing;
			globalDensityProfileBin = floor(globalDensityProfileBinInterval);
			globalDensityProfileBins[globalDensityProfileBin] += 1.0 / globalDensityProfileBinVolumes[globalDensityProfileBin];
		}

		//kinetic energy
		for (int a = 0; a < DIM; a++)
		{
			opd.vecKineticSumR1[a] = 0.0;
			opd.vecKineticSumI1[a] = 0.0;
		}

		int nCFD = 0;
		int paramOffset = 26;
		int offset = 0;
		for (auto& cfd : this->corrFuncData)
		{
			offset = nCFD * paramOffset;
			for (int a = 0; a < DIM; a++)
			{
				temp = (cfd.mcMillanSumD[n][a] + cfd.bcFactors[0][0] * cfd.splineSumsD[0][n][a] + cfd.bcFactors[0][1] * cfd.splineSumsD[1][n][a]);
				opd.vecKineticSumR1[a] += uR[0 + offset] * temp;
				opd.vecKineticSumI1[a] += uI[0 + offset] * temp;
				temp = (cfd.splineSumsD[2][n][a] + cfd.bcFactors[1][0] * cfd.splineSumsD[0][n][a] + cfd.bcFactors[1][1] * cfd.splineSumsD[1][n][a]);
				opd.vecKineticSumR1[a] += uR[1 + offset] * temp;
				opd.vecKineticSumI1[a] += uI[1 + offset] * temp;
			}
			temp = (cfd.mcMillanSumD2[n] + cfd.bcFactors[0][0] * cfd.splineSumsD2[0][n] + cfd.bcFactors[0][1] * cfd.splineSumsD2[1][n]);
			opd.kineticSumR2 += uR[0 + offset] * temp;
			opd.kineticSumI2 += uI[0 + offset] * temp;
			temp = (cfd.splineSumsD2[2][n] + cfd.bcFactors[1][0] * cfd.splineSumsD2[0][n] + cfd.bcFactors[1][1] * cfd.splineSumsD2[1][n]);
			opd.kineticSumR2 += uR[1 + offset] * temp;
			opd.kineticSumI2 += uI[1 + offset] * temp;

			for (int k = 2; k < paramOffset - 4; k++)
			{
				for (int a = 0; a < DIM; a++)
				{
					opd.vecKineticSumR1[a] += uR[k + offset] * cfd.splineSumsD[k + 1][n][a];
					opd.vecKineticSumI1[a] += uI[k + offset] * cfd.splineSumsD[k + 1][n][a];
				}
				opd.kineticSumR2 += uR[k + offset] * cfd.splineSumsD2[k + 1][n];
				opd.kineticSumI2 += uI[k + offset] * cfd.splineSumsD2[k + 1][n];
			}

			for (int a = 0; a < DIM; a++)
			{
				temp = (cfd.splineSumsD[cfd.numberOfSplines - 3][n][a] + cfd.bcFactors[2][0] * cfd.splineSumsD[cfd.numberOfSplines - 2][n][a] + cfd.bcFactors[2][1] * cfd.splineSumsD[cfd.numberOfSplines - 1][n][a]);
				opd.vecKineticSumR1[a] += uR[paramOffset - 4 + offset] * temp;
				opd.vecKineticSumI1[a] += uI[paramOffset - 4 + offset] * temp;
				temp = (cfd.constSumD[n][a] + cfd.bcFactors[3][0] * cfd.splineSumsD[cfd.numberOfSplines - 2][n][a] + cfd.bcFactors[3][1] * cfd.splineSumsD[cfd.numberOfSplines - 1][n][a]);
				opd.vecKineticSumR1[a] += uR[paramOffset - 3 + offset] * temp;
				opd.vecKineticSumI1[a] += uI[paramOffset - 3 + offset] * temp;
				temp = (cfd.linearSumD[n][a] + cfd.bcFactors[4][0] * cfd.splineSumsD[cfd.numberOfSplines - 2][n][a] + cfd.bcFactors[4][1] * cfd.splineSumsD[cfd.numberOfSplines - 1][n][a]);
				opd.vecKineticSumR1[a] += uR[paramOffset - 2 + offset] * temp;
				opd.vecKineticSumI1[a] += uI[paramOffset - 2 + offset] * temp;
			}
			temp = (cfd.splineSumsD2[cfd.numberOfSplines - 3][n] + cfd.bcFactors[2][0] * cfd.splineSumsD2[cfd.numberOfSplines - 2][n] + cfd.bcFactors[2][1] * cfd.splineSumsD2[cfd.numberOfSplines - 1][n]);
			opd.kineticSumR2 += uR[paramOffset - 4 + offset] * temp;
			opd.kineticSumI2 += uI[paramOffset - 4 + offset] * temp;
			temp = (cfd.constSumD2[n] + cfd.bcFactors[3][0] * cfd.splineSumsD2[cfd.numberOfSplines - 2][n] + cfd.bcFactors[3][1] * cfd.splineSumsD2[cfd.numberOfSplines - 1][n]);
			opd.kineticSumR2 += uR[paramOffset - 3 + offset] * temp;
			opd.kineticSumI2 += uI[paramOffset - 3 + offset] * temp;
			temp = (cfd.linearSumD2[n] + cfd.bcFactors[4][0] * cfd.splineSumsD2[cfd.numberOfSplines - 2][n] + cfd.bcFactors[4][1] * cfd.splineSumsD2[cfd.numberOfSplines - 1][n]);
			opd.kineticSumR2 += uR[paramOffset - 2 + offset] * temp;
			opd.kineticSumI2 += uI[paramOffset - 2 + offset] * temp;

			for (int a = 0; a < DIM; a++)
			{
				opd.vecKineticSumR1[a] += uR[paramOffset - 1 + offset] * cfd.logSumD[n][a];
				opd.vecKineticSumI1[a] += uI[paramOffset - 1 + offset] * cfd.logSumD[n][a];
			}
			opd.kineticSumR2 += uR[paramOffset - 1 + offset] * cfd.logSumD2[n];
			opd.kineticSumI2 += uI[paramOffset - 1 + offset] * cfd.logSumD2[n];

			nCFD++;
		}

		opd.kineticSumR1I1 = 2.0 * VectorDotProduct(opd.vecKineticSumR1, opd.vecKineticSumI1);
		opd.kineticSumR1 = VectorNorm2(opd.vecKineticSumR1);
		opd.kineticSumI1 = VectorNorm2(opd.vecKineticSumI1);

		kineticR += -pp.hbarOver2m * (opd.kineticSumR1 - opd.kineticSumI1 + opd.kineticSumR2);
		kineticI += -pp.hbarOver2m * (opd.kineticSumR1I1 + opd.kineticSumI2);

		kineticR_part1 += -pp.hbarOver2m * opd.kineticSumR1;
		kineticR_part2 += -pp.hbarOver2m * opd.kineticSumR2;
	}

	localEnergyR = kineticR + potentialIntern + potentialExtern;
	localEnergyI = kineticI;

	//cout << "potentialExtern=" << potentialExtern << endl;
	//cout << "potentialIntern=" << potentialIntern << endl;
	//cout << "kineticSumR1=" << kineticSumR1 << endl;
	//cout << "kineticSumI1=" << kineticSumI1 << endl;
	//cout << "kineticSumR2=" << kineticSumR2 << endl;
	//cout << "kineticSumI2=" << kineticSumI2 << endl;
	//cout << "kineticSumR1I1=" << kineticSumR1I1 << endl;

	int nCFD = 0;
	int paramOffset = 26;
	int offset = 0;
	for (auto& cfd : this->corrFuncData)
	{
		offset = nCFD * paramOffset;
		localOperators[0 + offset] = cfd.mcMillanSum + cfd.bcFactors[0][0] * cfd.splineSums[0] + cfd.bcFactors[0][1] * cfd.splineSums[1];
		localOperators[1 + offset] = cfd.splineSums[2] + cfd.bcFactors[1][0] * cfd.splineSums[0] + cfd.bcFactors[1][1] * cfd.splineSums[1];
		for (int i = 2; i < paramOffset - 4; i++)
		{
			localOperators[i + offset] = cfd.splineSums[i + 1];
		}
		localOperators[paramOffset - 4 + offset] = cfd.splineSums[cfd.numberOfSplines - 3] + cfd.bcFactors[2][0] * cfd.splineSums[cfd.numberOfSplines - 2] + cfd.bcFactors[2][1] * cfd.splineSums[cfd.numberOfSplines - 1];
		localOperators[paramOffset - 3 + offset] = cfd.constSum + cfd.bcFactors[3][0] * cfd.splineSums[cfd.numberOfSplines - 2] + cfd.bcFactors[3][1] * cfd.splineSums[cfd.numberOfSplines - 1];
		localOperators[paramOffset - 2 + offset] = cfd.linearSum + cfd.bcFactors[4][0] * cfd.splineSums[cfd.numberOfSplines - 2] + cfd.bcFactors[4][1] * cfd.splineSums[cfd.numberOfSplines - 1];
		localOperators[paramOffset - 1 + offset] = cfd.logSum;

		nCFD++;
	}

	for (int k = 0; k < N_PARAM; k++)
	{
		for (int j = 0; j < N_PARAM; j++)
		{
			localOperatorsMatrix[k][j] = localOperators[k] * localOperators[j];
		}
		localOperatorlocalEnergyR[k] = localOperators[k] * localEnergyR;
		localOperatorlocalEnergyI[k] = localOperators[k] * localEnergyI;
	}

	//otherExpectationValues[0] = kineticR;
	//otherExpectationValues[1] = potentialIntern;
	//otherExpectationValues[2] = exponent;
	otherExpectationValues[0] = kineticR_part1;
	otherExpectationValues[1] = kineticR_part2;
	otherExpectationValues[2] = kineticR;
	otherExpectationValues[3] = potentialIntern;
	otherExpectationValues[4] = wf;
	otherExpectationValues[5] = exponent;
	for (int i = 0; i < globalNumOfDensityProfileValues; i++)
	{
		//otherExpectationValues[3 + i] = globalDensityProfileBins[i];
	}
}

void BosonMixtureCluster::CalculateExpectationValues(ICorrelatedSamplingData* sample, vector<double>& uR, vector<double>& uI, double phiR, double phiI)
{
//TODO: implement
}

void BosonMixtureCluster::CalculateAdditionalSystemProperties(vector<vector<double> >& R, vector<double>& uR, vector<double>& uI, double phiR, double phiI)
{
	additionalObservables.ClearValues();
	vector<double> com;
	double rni;
	vector<double> vecrni(DIM);
	com = GetCenterOfMass(R);


	//r2
	double r;
	double r2Sum = 0.0;
	for (int i = 0; i < N; i++)
	{
		vector<double> diff = R[i] - com;
		r = VectorNorm(diff);
		r2Sum += r * r;
	}
	r2.value = r2Sum / (double) N;


	//angularDistribution
	double angle;

	//angle between particles 1-2-3
	angle = GetCornerAngle(R[0], R[1], R[2]);
	angularDistribution.AddToHistogram(0, angle, 1.0);

	//angle between particles 1-3-2
	angle = GetCornerAngle(R[0], R[2], R[1]);
	angularDistribution.AddToHistogram(1, angle, 1.0);

	//angle between particles 2-1-3
	angle = GetCornerAngle(R[1], R[0], R[2]);
	angularDistribution.AddToHistogram(2, angle, 1.0);


	//density from COM
	for (int i = 0; i < N; i++)
	{
		rni = VectorDisplacement(R[i], com, vecrni);
		if (rni < densityFromCOM.grid.max)
		{
			densityFromCOM.AddToHistogram(i, rni, 1.0);
		}
	}


	//distance matrix
	int index = 0;
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < i; j++)
		{
			rni = VectorDisplacement(R[i], R[j], vecrni);
			if (rni < particleDistances.grid.max)
			{
				particleDistances.AddToHistogram(index, rni, 1.0);
			}
			index++;
		}
	}
}

//void BosonMixtureCluster::CalculateAdditionalSystemProperties(vector<vector<double> >& R, vector<double>& uR, vector<double>& uI, double phiR, double phiI)
//{
//	ClearVector (additionalSystemProperties);
//	int index;
//
//	double r;
//	double r2Sum = 0.0;
//	vector<double> com;
//	com = GetCenterOfMass(R);
//	for (int i = 0; i < N; i++)
//	{
//		vector<double> diff = R[i] - com;
//		r = VectorNorm(diff);
//		r2Sum += r * r;
//	}
//	additionalSystemProperties[0] = r2Sum / (double) N;
//
//	double angle;
//	int anglebin;
//
//	//angle between particles 1-2-3
//	angle = GetCornerAngle(R[0], R[1], R[2]) / M_PI * 180.0;
//	anglebin = angle;
//	additionalSystemProperties[1 + anglebin] += 1;
//
//	//angle between particles 1-3-2
//	angle = GetCornerAngle(R[0], R[2], R[1]) / M_PI * 180.0;
//	anglebin = angle;
//	additionalSystemProperties[1 + 180 + anglebin] += 1;
//
//	//angle between particles 2-1-3
//	angle = GetCornerAngle(R[1], R[0], R[2]) / M_PI * 180.0;
//	anglebin = angle;
//	additionalSystemProperties[1 + 180 + 180 + anglebin] += 1;
//
//	double interval;
//	int bin;
//	double rni;
//	vector<double> vecrni(DIM);
//
//	//density from COM
//	vector<vector<double> > densityCOM;
//	InitVector(densityCOM, N, globalNumOfDensityProfileValues, 0.0);
//	for (int i = 0; i < N; i++)
//	{
//		rni = VectorDisplacement(R[i], com, vecrni);
//		if (rni < globalDensityProfileMaxDistance)
//		{
//			interval = rni / globalDensityProfileNodePointSpacing;
//			bin = floor(interval);
//			densityCOM[i][bin] += 1.0 / globalDensityProfileBinVolumes[bin];
//		}
//	}
//	index = 1 + 180 + 180 + 180;
//	for (int i = 0; i < N; i++)
//	{
//		for (int d = 0; d < globalNumOfDensityProfileValues; d++)
//		{
//			additionalSystemProperties[index] = densityCOM[i][d];
//			index++;
//		}
//	}
//
//	//distance matrix
//	vector<vector<vector<double> > > densityProfileMatrix;
//	InitVector(densityProfileMatrix, N, N, globalNumOfDensityProfileValues, 0.0);
//	for (int i = 0; i < N; i++)
//	{
//		for (int j = 0; j < i; j++)
//		{
//			rni = VectorDisplacement(R[i], R[j], vecrni);
//			if (rni < globalDensityProfileMaxDistance)
//			{
//				interval = rni / globalDensityProfileNodePointSpacing;
//				bin = floor(interval);
//				//densityProfileMatrix[i][j][bin] += 1.0 / globalDensityProfileBinVolumes[bin];
//				densityProfileMatrix[i][j][bin] += 1.0;
//			}
//		}
//	}
//	for (int i = 0; i < N; i++)
//	{
//		for (int j = 0; j < i; j++)
//		{
//			for (int d = 0; d < globalNumOfDensityProfileValues; d++)
//			{
//				additionalSystemProperties[index] = densityProfileMatrix[i][j][d];
//				index++;
//			}
//		}
//	}
//}

//void BosonMixtureCluster::GetCorrelationFunctionValues(vector<double>& uR, vector<double>& uI, double phiR, double phiI)
//{
//
//}

void BosonMixtureCluster::CalculateWavefunction(vector<vector<double> >& R, vector<double>& uR, vector<double>& uI, double phiR, double phiI)
{
	int ct;

	double sum = 0;
	int bin;
	double rni;
	double rni2;
	double rni3;
	vector<double> vecrni(DIM);

	for (auto& cfd : this->corrFuncData)
	{
		ClearVector(cfd.splineSums);
		cfd.mcMillanSum = 0;
		cfd.constSum = 0;
		cfd.logSum = 0;
		cfd.linearSum = 0;
	}

	for (int n = 0; n < N; n++)
	{
		for (int i = 0; i < n; i++)
		{
			ct = this->correlationTypes[n][i];
			CorrelationFunctionData& cfd = this->corrFuncData[ct];
			rni = VectorDisplacement(R[n], R[i], vecrni);
			if (rni <= cfd.rijSplit)
			{
				cfd.mcMillanSum += pow(rni, cfd.mcMillanFactor);
			}
			else if (rni >= cfd.rijTail)
			{
				cfd.constSum += 1.0;
				cfd.linearSum += rni;
			}
			else
			{
				auto interval = lower_bound(cfd.nodes.begin(), cfd.nodes.end(), rni);
				bin = interval - cfd.nodes.begin() - 1;
				rni2 = rni * rni;
				rni3 = rni2 * rni;

				for (int p = 0; p < 4; p++)
				{
					cfd.splineSums[bin - p] += cfd.splineWeights[bin - p][p][0] + cfd.splineWeights[bin - p][p][1] * rni + cfd.splineWeights[bin - p][p][2] * rni2 + cfd.splineWeights[bin - p][p][3] * rni3;
				}
			}
			cfd.logSum += log(rni);
		}
	}

	int nCFD = 0;
	int paramOffset = 26;
	int offset = 0;
	//for (auto& cfd : this->corrFuncData)
	//{
	//	WriteDataToFile(cfd.bcFactors, "bcFactors", "bcFactors");
	//}
	for (auto& cfd : this->corrFuncData)
	{
		offset = nCFD * paramOffset;
		//am
		sum += uR[0 + offset] * (cfd.mcMillanSum + cfd.bcFactors[0][0] * cfd.splineSums[0] + cfd.bcFactors[0][1] * cfd.splineSums[1]);
		//a2
		sum += uR[1 + offset] * (cfd.splineSums[2] + cfd.bcFactors[1][0] * cfd.splineSums[0] + cfd.bcFactors[1][1] * cfd.splineSums[1]);
		for (int i = 2; i < paramOffset - 4; i++)
		{
			sum += uR[i + offset] * cfd.splineSums[i + 1];
		}
		//a(n-2)
		sum += uR[paramOffset - 4 + offset] * (cfd.splineSums[cfd.numberOfSplines - 3] + cfd.bcFactors[2][0] * cfd.splineSums[cfd.numberOfSplines - 2] + cfd.bcFactors[2][1] * cfd.splineSums[cfd.numberOfSplines - 1]);
		//ac
		sum += uR[paramOffset - 3 + offset] * (cfd.constSum + cfd.bcFactors[3][0] * cfd.splineSums[cfd.numberOfSplines - 2] + cfd.bcFactors[3][1] * cfd.splineSums[cfd.numberOfSplines - 1]);
		//alin
		sum += uR[paramOffset - 2 + offset] * (cfd.linearSum + cfd.bcFactors[4][0] * cfd.splineSums[cfd.numberOfSplines - 2] + cfd.bcFactors[4][1] * cfd.splineSums[cfd.numberOfSplines - 1]);

		sum += uR[paramOffset - 1 + offset] * cfd.logSum;
		nCFD++;
	}

	wf = exp(sum + phiR);
	exponent = sum;

//cout << "sum=" << sum << "\t\tsumTmp=" << sumTmp << "\t\twf=" << value << endl;
//cout << "sum=" << sum << "\t\twf=" << wf << "\t\tlogSum=" << logSum << "\t\tlinearSum=" << linearSum << "\t\tphiR=" << phiR << endl;
}

void BosonMixtureCluster::CalculateWavefunction(ICorrelatedSamplingData* sample, vector<double>& uR, vector<double>& uI, double phiR, double phiI)
{
//TODO: implement
}

void BosonMixtureCluster::CalculateWFChange(vector<vector<double> >& R, vector<double>& uR, vector<double>& uI, double phiR, double phiI, int changedParticleIndex, vector<double>& oldPosition)
{
	int ct;

	double sum = 0;
	int bin;
	double rni;
	double rni2;
	double rni3;
	vector<double> vecrni(DIM);

	for (auto& cfd : this->corrFuncData)
	{
		ClearVector(cfd.sumOldPerBin);
		ClearVector(cfd.sumNewPerBin);
		cfd.mcMillanOld = 0;
		cfd.mcMillanNew = 0;
		cfd.constOld = 0;
		cfd.constNew = 0;
		cfd.logOld = 0;
		cfd.logNew = 0;
		cfd.linearOld = 0;
		cfd.linearNew = 0;
	}

	for (int i = 0; i < N; i++)
	{
		if (i != changedParticleIndex)
		{
			ct = this->correlationTypes[i][changedParticleIndex];
			CorrelationFunctionData& cfd = this->corrFuncData[ct];
			rni = VectorDisplacement(R[i], oldPosition, vecrni);
			if (rni < cfd.rijSplit)
			{
				cfd.mcMillanOld += pow(rni, cfd.mcMillanFactor);
			}
			else if (rni >= cfd.rijTail)
			{
				cfd.constOld += 1.0;
				cfd.linearOld += rni;
			}
			else
			{
				auto interval = lower_bound(cfd.nodes.begin(), cfd.nodes.end(), rni);
				bin = interval - cfd.nodes.begin() - 1;
				rni2 = rni * rni;
				rni3 = rni2 * rni;

				for (int p = 0; p < 4; p++)
				{
					cfd.sumOldPerBin[bin - p] += cfd.splineWeights[bin - p][p][0] + cfd.splineWeights[bin - p][p][1] * rni + cfd.splineWeights[bin - p][p][2] * rni2 + cfd.splineWeights[bin - p][p][3] * rni3;
				}
			}
			cfd.logOld += log(rni);

			rni = VectorDisplacement(R[i], R[changedParticleIndex], vecrni);
			if (rni < cfd.rijSplit)
			{
				cfd.mcMillanNew += pow(rni, cfd.mcMillanFactor);
			}
			else if (rni >= cfd.rijTail)
			{
				cfd.constNew += 1.0;
				cfd.linearNew += rni;
			}
			else
			{
				auto interval = lower_bound(cfd.nodes.begin(), cfd.nodes.end(), rni);
				bin = interval - cfd.nodes.begin() - 1;
				rni2 = rni * rni;
				rni3 = rni2 * rni;

				for (int p = 0; p < 4; p++)
				{
					cfd.sumNewPerBin[bin - p] += cfd.splineWeights[bin - p][p][0] + cfd.splineWeights[bin - p][p][1] * rni + cfd.splineWeights[bin - p][p][2] * rni2 + cfd.splineWeights[bin - p][p][3] * rni3;
				}
			}
			cfd.logNew += log(rni);
		}
	}

	for (auto& cfd : this->corrFuncData)
	{
		cfd.mcMillanSumNew = max(0.0, cfd.mcMillanSum - cfd.mcMillanOld + cfd.mcMillanNew);
		for (int i = 0; i < cfd.numberOfSplines; i++)
		{
			cfd.splineSumsNew[i] = max(0.0, cfd.splineSums[i] - cfd.sumOldPerBin[i] + cfd.sumNewPerBin[i]);
		}
		cfd.constSumNew = cfd.constSum - cfd.constOld + cfd.constNew;
		cfd.logSumNew = cfd.logSum - cfd.logOld + cfd.logNew;
		cfd.linearSumNew = max(0.0, cfd.linearSum - cfd.linearOld + cfd.linearNew);
	}
	int nCFD = 0;
	int paramOffset = 26;
	int offset = 0;
	for (auto& cfd : this->corrFuncData)
	{
		offset = nCFD * paramOffset;
		//am
		sum += uR[0 + offset] * (cfd.mcMillanSumNew + cfd.bcFactors[0][0] * cfd.splineSumsNew[0] + cfd.bcFactors[0][1] * cfd.splineSumsNew[1]);
		//a2
		sum += uR[1 + offset] * (cfd.splineSumsNew[2] + cfd.bcFactors[1][0] * cfd.splineSumsNew[0] + cfd.bcFactors[1][1] * cfd.splineSumsNew[1]);
		for (int i = 2; i < paramOffset - 4; i++)
		{
			sum += uR[i + offset] * cfd.splineSumsNew[i + 1];
		}
		//a(n-2)
		sum += uR[paramOffset - 4 + offset] * (cfd.splineSumsNew[cfd.numberOfSplines - 3] + cfd.bcFactors[2][0] * cfd.splineSumsNew[cfd.numberOfSplines - 2] + cfd.bcFactors[2][1] * cfd.splineSumsNew[cfd.numberOfSplines - 1]);
		//ac
		sum += uR[paramOffset - 3 + offset] * (cfd.constSumNew + cfd.bcFactors[3][0] * cfd.splineSumsNew[cfd.numberOfSplines - 2] + cfd.bcFactors[3][1] * cfd.splineSumsNew[cfd.numberOfSplines - 1]);
		//alin
		sum += uR[paramOffset - 2 + offset] * (cfd.linearSumNew + cfd.bcFactors[4][0] * cfd.splineSumsNew[cfd.numberOfSplines - 2] + cfd.bcFactors[4][1] * cfd.splineSumsNew[cfd.numberOfSplines - 1]);

		sum += uR[paramOffset - 1 + offset] * cfd.logSumNew;
		nCFD++;
	}
	wfNew = exp(sum + phiR);
	exponentNew = sum;
}

double BosonMixtureCluster::CalculateWFQuotient(vector<vector<double> >& R, vector<double>& uR, vector<double>& uI, double phiR, double phiI, int changedParticleIndex, vector<double>& oldPosition)
{
	CalculateWFChange(R, uR, uI, phiR, phiI, changedParticleIndex, oldPosition);
	double wfQuotient = exp(2.0 * (exponentNew - exponent));

//cout << "qu=" << wfQuotient << endl;

	return wfQuotient;
}

void BosonMixtureCluster::AcceptMove()
{
	wf = wfNew;
	exponent = exponentNew;
	for (auto& cfd : this->corrFuncData)
	{
		cfd.mcMillanSum = cfd.mcMillanSumNew;
		cfd.splineSums = cfd.splineSumsNew;
		cfd.constSum = cfd.constSumNew;
		cfd.logSum = cfd.logSumNew;
		cfd.linearSum = cfd.linearSumNew;
	}
}

void BosonMixtureCluster::InitCorrelatedSamplingData(vector<ICorrelatedSamplingData*>& data)
{
}

void BosonMixtureCluster::FillCorrelatedSamplingData(ICorrelatedSamplingData* data)
{
//TODO: implement
}

}
