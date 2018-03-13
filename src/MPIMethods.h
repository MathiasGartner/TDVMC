#pragma once

#include "mpi.h"
#include <vector>

using namespace std;

namespace MPIMethods
{

int numOfProcesses;
int processRank;
int rootRank;

void BroadcastValue(double* data)
{
	MPI_Bcast(data, 1, MPI_DOUBLE, rootRank, MPI_COMM_WORLD);
}

void BroadcastValues(double* data, int count)
{
	MPI_Bcast(data, count, MPI_DOUBLE, rootRank, MPI_COMM_WORLD);
}

void BroadcastValues(vector<double>& data)
{
	double* valueArray = new double[data.size()];
	for (unsigned int i = 0; i < data.size(); i++)
	{
		valueArray[i] = data[i];
	}
	MPI_Bcast(valueArray, data.size(), MPI_DOUBLE, rootRank, MPI_COMM_WORLD);
	for (unsigned int i = 0; i < data.size(); i++)
	{
		data[i] = valueArray[i];
	}
	delete[] valueArray;
}

void BroadcastValues(vector<vector<double> >& data)
{
	for (unsigned int i = 0; i < data.size(); i++)
	{
		BroadcastValues(data[i]);
	}
}

vector<double> ReduceToMinMaxMean(double data)
{
	vector<double> result = { };
	double own = 0;
	double min = 0;
	double max = 0;
	double sum = 0;

	own = data;
	MPI_Reduce(&own, &min, 1, MPI_DOUBLE, MPI_MIN, rootRank, MPI_COMM_WORLD);
	MPI_Reduce(&own, &max, 1, MPI_DOUBLE, MPI_MAX, rootRank, MPI_COMM_WORLD);
	MPI_Reduce(&own, &sum, 1, MPI_DOUBLE, MPI_SUM, rootRank, MPI_COMM_WORLD);
	if (processRank == rootRank)
	{
		result = { min, max, sum / (double)numOfProcesses };
	}
	return result;
}

void ReduceToAverage(double* data)
{
	double ownValues;
	double reducedValues;

	ownValues = *data;
	MPI_Reduce(&ownValues, &reducedValues, 1, MPI_DOUBLE, MPI_SUM, rootRank, MPI_COMM_WORLD);
	if (processRank == rootRank)
	{
		*data = reducedValues / (double) numOfProcesses;
	}
}

void ReduceToAverage(vector<double>& data)
{
	double* ownValues = new double[data.size()];
	double* reducedValues = new double[data.size()];

	for (unsigned int i = 0; i < data.size(); i++)
	{
		ownValues[i] = data[i];
	}
	MPI_Reduce(ownValues, reducedValues, data.size(), MPI_DOUBLE, MPI_SUM, rootRank, MPI_COMM_WORLD);
	//INFO: only root process holds the average values from all processes.
	if (processRank == rootRank)
	{
		for (unsigned int i = 0; i < data.size(); i++)
		{
			data[i] = reducedValues[i] / (double) numOfProcesses;
		}
	}
	delete[] ownValues;
	delete[] reducedValues;
}

void ReduceToAverage(vector<vector<double> >& data)
{
	for (unsigned int i = 0; i < data.size(); i++)
	{
		ReduceToAverage(data[i]);
	}
}

}
