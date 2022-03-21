#include "SplinedFunction.h"

namespace WFParts
{

SplinedFunction::SplinedFunction(int splineOrder_)
{
	splineOrder = splineOrder_;
	numOfSplineParts = splineOrder + 1;

	nodeSpacing = 0.0;

	numberOfSplines = 0;

	numberOfSpecialParametersStart = 0;
	numberOfSpecialParametersEnd = 0;
	numberOfStandardParameters = 0;
	numberOfTotalParameters = 0;
	np1 = 0;
	np2 = 0;
	np3 = 0;
}

void SplinedFunction::Init()
{
	InitVector(splineSums, numberOfSplines, 0.0);
	InitVector(splineSumsD, numberOfSplines, N, DIM, 0.0);
	InitVector(splineSumsD2, numberOfSplines, N, 0.0);

	InitVector(splineSumsNew, numberOfSplines, 0.0);
	InitVector(sumOldPerBin, numberOfSplines, 0.0);
	InitVector(sumNewPerBin, numberOfSplines, 0.0);
}

int SplinedFunction::BinIndex(double r)
{
	return GetBinIndex(this->nodes, r);
}

}
