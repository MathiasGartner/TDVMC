/*-----------------------------------------------------------------------------
 *
 * 		Name:			CorrelationFunctionData.cpp
 * 		Author:			Mathias Gartner
 * 		Description:	see CorrelationFunctionData.h
 *
 * --------------------------------------------------------------------------*/

#include "CorrelationFunctionData.h"

CorrelationFunctionData::CorrelationFunctionData()
{
	numberOfSplines = 0;

	mcMillanSum = 0;
	constSum = 0;
	logSum = 0;
	linearSum = 0;

	rijSplit = 0;
	rijTail = 0;
	mcMillanFactor = 0.0;

	mcMillanSumNew = 0;
	constSumNew = 0;
	logSumNew = 0;
	linearSumNew = 0;

	mcMillanOld = 0;
	mcMillanNew = 0;
	constOld = 0;
	constNew = 0;
	logOld = 0;
	logNew = 0;
	linearOld = 0;
	linearNew = 0;

}

void CorrelationFunctionData::Init()
{
	mcMillanSum = 0;
	mcMillanSumD.resize(N);
	for (auto &n : mcMillanSumD)
	{
		n.resize(DIM);
	}
	ClearVector(mcMillanSumD);
	mcMillanSumD2.resize(N);
	ClearVector(mcMillanSumD2);
	constSum = 0;
	constSumD.resize(N);
	for (auto &n : constSumD)
	{
		n.resize(DIM);
	}
	ClearVector(constSumD);
	constSumD2.resize(N);
	ClearVector(constSumD2);
	logSum = 0;
	logSumD.resize(N);
	for (auto &n : logSumD)
	{
		n.resize(DIM);
	}
	ClearVector(logSumD);
	logSumD2.resize(N);
	ClearVector(logSumD2);
	linearSum = 0;
	linearSumD.resize(N);
	for (auto &n : linearSumD)
	{
		n.resize(DIM);
	}
	ClearVector(linearSumD);
	linearSumD2.resize(N);
	ClearVector(linearSumD2);
	splineSums.resize(numberOfSplines);
	splineSumsD.resize(numberOfSplines);
	for (auto &k : splineSumsD)
	{
		k.resize(N);
		for (auto &n : k)
		{
			n.resize(DIM);
		}
	}
	splineSumsD2.resize(numberOfSplines);
	for (auto &k : splineSumsD2)
	{
		k.resize(N);
	}

	mcMillanSumNew = 0.0;
	constSumNew = 0.0;
	logSumNew = 0.0;
	linearSumNew = 0.0;
	splineSumsNew.resize(numberOfSplines);
	sumOldPerBin.resize(numberOfSplines);
	sumNewPerBin.resize(numberOfSplines);
}
