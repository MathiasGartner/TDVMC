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
	np1 = 0;
	np2 = 0;
	np3 = 0;

}

}
