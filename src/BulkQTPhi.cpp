#include "BulkQTPhi.h"

BulkQTPhi::BulkQTPhi(string configDirectory) : BulkQT(configDirectory)
{
	this->USE_NORMALIZATION_AND_PHASE = true;
}

void BulkQTPhi::CalculateWavefunction(vector<vector<double> >& R, vector<double>& uR, vector<double>& uI, double phiR, double phiI)
{
	BulkQT::CalculateWavefunction(R, uR, uI, phiR, phiI);
	wf *= exp(phiR);
}

void BulkQTPhi::CalculateWFChange(vector<vector<double> >& R, vector<double>& uR, vector<double>& uI, double phiR, double phiI, int changedParticleIndex, vector<double>& oldPosition)
{
	BulkQT::CalculateWFChange(R, uR, uI, phiR, phiI, changedParticleIndex, oldPosition);
	wfNew *= exp(phiR);
}
