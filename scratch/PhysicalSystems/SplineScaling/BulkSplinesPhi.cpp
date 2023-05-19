#include "BulkSplinesPhi.h"

namespace PhysicalSystems
{

BulkSplinesPhi::BulkSplinesPhi(vector<double>& params, string configDirectory) :
		BulkSplines(params, configDirectory)
{
	this->USE_NORMALIZATION_AND_PHASE = true;
}

void BulkSplinesPhi::CalculateWavefunction(vector<vector<double> >& R, vector<double>& uR, vector<double>& uI, double phiR, double phiI)
{
	BulkSplines::CalculateWavefunction(R, uR, uI, phiR, phiI);
	wf = exp(exponent + phiR);
}

void BulkSplinesPhi::CalculateWFChange(vector<vector<double> >& R, vector<double>& uR, vector<double>& uI, double phiR, double phiI, int changedParticleIndex, vector<double>& oldPosition)
{
	BulkSplines::CalculateWFChange(R, uR, uI, phiR, phiI, changedParticleIndex, oldPosition);
	wfNew = exp(exponentNew + phiR);
}

}
