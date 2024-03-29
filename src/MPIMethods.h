/*-----------------------------------------------------------------------------
 *
 * 		Name:			MPIMethods.h
 * 		Author:			Mathias Gartner
 * 		Description:	Wrapper functions for the standard MPI methods as well
 * 						as functions combining multiple MPI calls.
 *
 * --------------------------------------------------------------------------*/

#pragma once

#include "mpi.h"

#include "Observables/Observable.h"
#include "Observables/ObservableCollection.h"
#include "Observables/ObservableV.h"
#include "Observables/ObservableVsOnGrid.h"

#include <cstring>
#include <vector>

using namespace std;

namespace MPIMethods
{

int numOfProcesses;
int processRank;
int rootRank;
bool isRootRank;

void Barrier()
{
	MPI_Barrier(MPI_COMM_WORLD);
}

void BroadcastValue(string* data, int maxLength)
{
	char tmp[maxLength];
	strcpy(tmp, (*data).c_str());
	MPI_Bcast(tmp, maxLength, MPI_CHAR, rootRank, MPI_COMM_WORLD);
	*data = string(tmp);
}

void BroadcastValue(int* data)
{
	MPI_Bcast(data, 1, MPI_INT, rootRank, MPI_COMM_WORLD);
}

void BroadcastValues(vector<int>& data)
{
	MPI_Bcast(data.data(), data.size(), MPI_INT, rootRank, MPI_COMM_WORLD);
}

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
	MPI_Bcast(data.data(), data.size(), MPI_DOUBLE, rootRank, MPI_COMM_WORLD);
}

void BroadcastValues(vector<vector<double> >& data)
{
	//TODO: if needed in TDVMC-calculations -> check performance
	for (unsigned int i = 0; i < data.size(); i++)
	{
		BroadcastValues(data[i]);
	}
}

void BroadcastValues(vector<vector<vector<double> > >& data)
{
	//TODO: if needed in TDVMC-calculations -> check performance
	for (unsigned int i = 0; i < data.size(); i++)
	{
		for (unsigned int j = 0; j < data.size(); j++)
		{
			BroadcastValues(data[i][j]);
		}
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
	if (isRootRank)
	{
		result =
		{	min, max, sum / (double)numOfProcesses};
	}
	return result;
}

bool IsAnyTrue(bool value)
{
	int ownValue = value ? 1 : 0;
	int reducedValue;
	MPI_Allreduce(&ownValue, &reducedValue, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
	return reducedValue > 0;
}

bool IsAnyFalse(bool value)
{
	int ownValue = value ? 1 : 0;
	int reducedValue;
	MPI_Allreduce(&ownValue, &reducedValue, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
	return reducedValue < 1;
}

vector<bool> CollectValues(bool value)
{
	vector<bool> values(numOfProcesses);
	int ownValue = value ? 1 : 0;
	int* gatheredValues = new int[numOfProcesses];
	MPI_Gather(&ownValue, 1, MPI_INT, gatheredValues, 1, MPI_INT, rootRank, MPI_COMM_WORLD);
	for (int i = 0; i < numOfProcesses; i++)
	{
		values[i] = gatheredValues[i] == 1;
	}
	delete[] gatheredValues;
	return values;
}

void ReduceToAverage(double* data)
{
	double ownValues;
	double reducedValues;

	ownValues = *data;
	MPI_Reduce(&ownValues, &reducedValues, 1, MPI_DOUBLE, MPI_SUM, rootRank, MPI_COMM_WORLD);
	if (isRootRank)
	{
		*data = reducedValues / (double) numOfProcesses;
	}
}

double ReduceToAverage(int* data)
{
	double average = 0;
	int ownValues;
	int reducedValues;

	ownValues = *data;
	MPI_Reduce(&ownValues, &reducedValues, 1, MPI_INT, MPI_SUM, rootRank, MPI_COMM_WORLD);
	if (isRootRank)
	{
		average = reducedValues / (double) numOfProcesses;
		*data = average;
	}
	return average;
}

double ReduceToAverage(long long* data)
{
	double average = 0;
	long long ownValues;
	long long reducedValues;

	ownValues = *data;
	MPI_Reduce(&ownValues, &reducedValues, 1, MPI_LONG_LONG_INT, MPI_SUM, rootRank, MPI_COMM_WORLD);
	if (isRootRank)
	{
		average = reducedValues / (double) numOfProcesses;
		*data = average;
	}
	return average;
}

void ReduceToAverage(double* data, int count)
{
	double* reducedValues = new double[count];
	MPI_Reduce(data, reducedValues, count, MPI_DOUBLE, MPI_SUM, rootRank, MPI_COMM_WORLD);
	if (isRootRank)
	{
		for (int i = 0; i < count; i++)
		{
			data[i] = reducedValues[i] / (double) numOfProcesses;
		}
	}
	delete[] reducedValues;
}

void ReduceToAverage(vector<double>& data)
{
	double* reducedValues = new double[data.size()]; //TODO: check on zusie what method is faster

	MPI_Reduce(data.data(), reducedValues, data.size(), MPI_DOUBLE, MPI_SUM, rootRank, MPI_COMM_WORLD);
	//INFO: only root process holds the average values from all processes.
	if (isRootRank)
	{
		for (unsigned int i = 0; i < data.size(); i++)
		{
			data[i] = reducedValues[i] / (double) numOfProcesses;
		}
	}

	delete[] reducedValues;
}

void ReduceToAverage(Observables::Observable& data)
{
	double ownValues;
	double reducedValues;

	ownValues = data.value;
	MPI_Reduce(&ownValues, &reducedValues, 1, MPI_DOUBLE, MPI_SUM, rootRank, MPI_COMM_WORLD);
	if (isRootRank)
	{
		data.value = reducedValues / (double) numOfProcesses;
	}
}

void ReduceToAverage(Observables::ObservableV& data)
{
	double* reducedValues = new double[data.values.size()]; //TODO: check on zusie what method is faster

	MPI_Reduce(data.values.data(), reducedValues, data.values.size(), MPI_DOUBLE, MPI_SUM, rootRank, MPI_COMM_WORLD);
	//INFO: only root process holds the average values from all processes.
	if (isRootRank)
	{
		for (unsigned int i = 0; i < data.values.size(); i++)
		{
			data.values[i] = reducedValues[i] / (double) numOfProcesses;
		}
	}

	delete[] reducedValues;
}

void ReduceToAverage(Observables::ObservableVsOnGrid& data)
{
	if (data.observablesV.size() > 0)
	{
		double* reducedValues = new double[data.grid.count]; //TODO: check on zusie what method is faster

		for (unsigned int j = 0; j < data.observablesV.size(); j++)
		{
			//TODO: just use ReduceToAverage(data.observablesV[j]); ??
			MPI_Reduce(data.observablesV[j].values.data(), reducedValues, data.observablesV[j].values.size(), MPI_DOUBLE, MPI_SUM, rootRank, MPI_COMM_WORLD);
			//INFO: only root process holds the average values from all processes.
			if (isRootRank)
			{
				for (int i = 0; i < data.grid.count; i++)
				{
					data.observablesV[j].values[i] = reducedValues[i] / (double) numOfProcesses;
				}
			}
		}

		delete[] reducedValues;
	}
}

void ReduceToAverage(Observables::ObservableVsOnMultiGrid& data)
{
	if (data.observablesV.size() > 0)
	{
		double* reducedValues = new double[data.totalGridPoints]; //TODO: check on zusie what method is faster

		for (unsigned int j = 0; j < data.observablesV.size(); j++)
		{
			//TODO: just use ReduceToAverage(data.observablesV[j]); ??
			MPI_Reduce(data.observablesV[j].values.data(), reducedValues, data.observablesV[j].values.size(), MPI_DOUBLE, MPI_SUM, rootRank, MPI_COMM_WORLD);
			//INFO: only root process holds the average values from all processes.
			if (isRootRank)
			{
				for (int i = 0; i < data.totalGridPoints; i++)
				{
					data.observablesV[j].values[i] = reducedValues[i] / (double) numOfProcesses;
				}
			}
		}

		delete[] reducedValues;
	}
}

void ReduceToAverage(Observables::ObservableCollection& data)
{
	for (auto& obs : data.observables)
	{
		if (auto o = dynamic_cast<Observables::ObservableVsOnMultiGrid*>(obs))
		{
			ReduceToAverage(*o);
		}
		else if (auto o = dynamic_cast<Observables::ObservableVsOnGrid*>(obs))
		{
			ReduceToAverage(*o);
		}
		else if (auto o = dynamic_cast<Observables::ObservableV*>(obs))
		{
			ReduceToAverage(*o);
		}
		else if (auto o = dynamic_cast<Observables::Observable*>(obs))
		{
			ReduceToAverage(*o);
		}
		else
		{
			throw std::invalid_argument("ReduceToAverage of IObservable* not implemented for given type.");
		}
	}
}

//INFO: is a little bit slower than ReduceToAverage with double* reducedValues
//void ReduceToAverage(vector<double>& data)
//{
//	vector<double> reducedValues(data.size());
//
//	MPI_Reduce(data.data(), reducedValues.data(), data.size(), MPI_DOUBLE, MPI_SUM, rootRank, MPI_COMM_WORLD);
//	//INFO: only root process holds the average values from all processes.
//	if (isRootRank)
//	{
//		for (unsigned int i = 0; i < data.size(); i++)
//		{
//			data[i] = reducedValues[i] / (double) numOfProcesses;
//		}
//	}
//}

void ReduceToAverage(vector<vector<double> >& data)
{
	int size1 = data.size();
	int size2 = size1 > 0 ? data[0].size() : 0;
	int totalSize = size1 * size2;
	//TODO: besseres kriterium überlegen und auch auf anderen servern testen
	//		der wert 500 ist für zusie bie 256 prozessen.
	//		wie viele elemente insgesamt? quadratische matrix oder eher recht lange/breite matrix? ...
	if (size2 > 500)
	{
		//TODO: warum funtioniert diese methode nicht? ist extrem langsam wenn auf zusie 257 prozesse gestartet werden (für 256 noch ganz normal)
		for (unsigned int i = 0; i < data.size(); i++)
		{
			ReduceToAverage(data[i]);
		}
	}
	else
	{
		//TODO: funktioniert auch für >257 prozesse, ist aber langsamer. effizientere methode möglich?
		vector<double> combinedValues(totalSize);
		for (int i = 0; i < totalSize; i++)
		{
			combinedValues[i] = data[i / size2][i % size2];
		}
		ReduceToAverage(combinedValues);
		if (isRootRank)
		{
			for (int i = 0; i < totalSize; i++)
			{
				data[i / size2][i % size2] = combinedValues[i];
			}
		}
	}
}

map<int, int> GatherHistogram(int data)
{
	vector<int> gatheredValues(numOfProcesses);
	map<int, int> histogram;

	MPI_Gather(&data, 1, MPI_INT, gatheredValues.data(), 1, MPI_INT, rootRank, MPI_COMM_WORLD);
	if (isRootRank)
	{
		for (int i = 0; i < numOfProcesses; i++)
		{
			histogram[gatheredValues[i]]++;
		}
	}
	return histogram;
}

void GetCPUAllocation(bool printOnlyMultipleAllocations)
{
	bool print = false;
	bool infoPrinted = false;
	int id = get_cpu_id();
	map<int, int> histogram = MPIMethods::GatherHistogram(id);

	for (auto const& h : histogram)
	{
		if (h.second > 1)
		{
			if (printOnlyMultipleAllocations)
			{
				if (!infoPrinted)
				{
					cout << "########################" <<  endl;
					cout << "muliple CPU allocations:" <<  endl;
					infoPrinted = true;
				}
				cout << h.first << ":" << h.second << endl;
			}
			else
			{
				print = true;
				break;
			}
		}
	}
	if (infoPrinted)
	{
		cout << "########################" <<  endl;
	}
	if (print)
	{
		for (auto const& h : histogram)
		{
			cout << h.first << ":" << h.second << endl;
		}
	}
}

}
