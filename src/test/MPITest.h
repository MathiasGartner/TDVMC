#pragma once

#include "mpi.h"

#include "../Constants.h"
#include "../MPIMethods.h"
#include "../Timer.h"
#include "../Utils.h"

using namespace std;

namespace Test
{

int numOfProcesses = 1;
int rootRank = 0;
int processRank = 0;

void CreateRandomArray(double* arr, int n)
{
	mt19937_64 generator(time(NULL));
	uniform_real_distribution<double> distUniform(0.0, 1.0);
	for (int i = 0; i < n; i++)
	{
		arr[i] = distUniform(generator);
	}
}

void CreateRandomVector(vector<double>& v, int n, int randType)
{
	v.resize(n);
	if (randType == 0)
	{
		for (int i = 0; i < n; i++)
		{
			v[i] = rand() / 1.0;
		}
	}
	if (randType == 1)
	{
		mt19937 generator(time(NULL));
		uniform_real_distribution<double> distUniform(0.0, 1.0);
		for (int i = 0; i < n; i++)
		{
			v[i] = distUniform(generator);
		}
	} else if (randType == 2)
	{
		mt19937_64 generator(time(NULL));
		uniform_real_distribution<double> distUniform(0.0, 1.0);
		for (int i = 0; i < n; i++)
		{
			v[i] = distUniform(generator);
		}
	}
}

void CreateRandomVector2D(vector<vector<double> >& v, int n, int randType)
{
	v.resize(n);
	for (int i = 0; i < n; i++)
	{
		CreateRandomVector(v[i], n, randType);
	}
}

void PrintTimings(vector<double> timings, string message)
{
	if (processRank == rootRank)
	{
		cout << message << endl;
		cout << "duration: min = " << to_string(timings[0]) << " ms" << endl;
		cout << "          max = " << to_string(timings[1]) << " ms" << endl;
		cout << "          <t> = " << to_string(timings[2]) << " ms" << endl;
	}
}

void GetCPUUsage()
{
	bool print = false;
	int id = get_cpu_id();
	map<int, int> histogram = MPIMethods::GatherHistogram(id);

	for (auto const& h : histogram)
	{
		if (h.second > 1)
		{
			print = true;
		}
	}
	if (print)
	{
		for (auto const& h : histogram)
		{
			cout << h.first << ":" << h.second << endl;
		}
	}
}

void TestMPIData()
{
	string str;
	int intValue;
	double value;
	double* data = new double[3];
	vector<double> v;
	vector<vector<double> > v2;
	map<int, int> h;

	//Broadcast string test
	str = "test";
	if (processRank == rootRank)
	{
		str = "broadcast!";
	}
	MPIMethods::BroadcastValue(&str, 300);
	if (str != "broadcast!")
	{
		cout << "ERROR Broadcast string on process " << to_string(processRank) << "!!" << endl;
	}
	else if (processRank == rootRank)
	{
		cout << "Broadcast string done" << endl;
	}

	//Broadcast int test
	intValue = 1;
	if (processRank == rootRank)
	{
		intValue = 3;
	}
	MPIMethods::BroadcastValue(&intValue);
	if (intValue != 3)
	{
		cout << "ERROR Broadcast int on process " << to_string(processRank) << "!!" << endl;
	}
	else if (processRank == rootRank)
	{
		cout << "Broadcast int done" << endl;
	}

	//Broadcast double test
	value = 1.0;
	if (processRank == rootRank)
	{
		value = 3.3;
	}
	MPIMethods::BroadcastValue(&value);
	if (value != 3.3)
	{
		cout << "ERROR on process " << to_string(processRank) << "!!" << endl;
	}
	else if (processRank == rootRank)
	{
		cout << "Broadcast double done" << endl;
	}

	//Broadcast double array test
	data[0] = 1.0; data[1] = 2.0; data[2] = 3.0;
	if (processRank == rootRank)
	{
		data[0] = 1.1; data[1] = 2.2; data[2] = 3.3;
	}
	MPIMethods::BroadcastValues(data, 3);
	if (data[0] != 1.1 ||
		data[1] != 2.2 ||
		data[2] != 3.3)
	{
		cout << "ERROR on process " << to_string(processRank) << "!!" << endl;
	}
	else if (processRank == rootRank)
	{
		cout << "Broadcast double array done" << endl;
	}

	//Broadcast double vector test
	v = {1.0, 2.0, 3.0};
	if (processRank == rootRank)
	{
		v = {1.1, 2.2, 3.3};
	}
	MPIMethods::BroadcastValues(v);
	if (v[0] != 1.1 ||
		v[1] != 2.2 ||
		v[2] != 3.3)
	{
		cout << "ERROR on process " << to_string(processRank) << "!!" << endl;
	}
	else if (processRank == rootRank)
	{
		cout << "Broadcast double vector done" << endl;
	}

	//Broadcast double 2D vector test
	v2 = {{1.0, 2.0, 3.0}, {4.0, 5.0, 6.0}};
	if (processRank == rootRank)
	{
		v2 = {{1.1, 2.2, 3.3}, {4.4, 5.5, 6.6}};
	}
	MPIMethods::BroadcastValues(v2);
	if (v2[0][0] != 1.1 ||
		v2[0][1] != 2.2 ||
		v2[0][2] != 3.3 ||
		v2[1][0] != 4.4 ||
		v2[1][1] != 5.5 ||
		v2[1][2] != 6.6)
	{
		cout << "ERROR on process " << to_string(processRank) << "!!" << endl;
	}
	else if (processRank == rootRank)
	{
		cout << "Broadcast double 2D vector done" << endl;
	}

	//ReduceToMinMaxMean test
	value = processRank * 2.0;
	v = MPIMethods::ReduceToMinMaxMean(value);
	if (processRank == rootRank)
	{
		if (v[0] != 0.0 ||
			v[1] != 2.0 * (numOfProcesses - 1) ||
			v[2] != (numOfProcesses - 1))
		{
			cout << "ERROR calculating ReduceToMinMaxMean!!" << endl;
		}
		cout << "ReduceToMinMaxMean done" << endl;
	}

	//ReduceToAverage double test
	value = processRank * 2.0;
	MPIMethods::ReduceToAverage(&value);
	if (processRank == rootRank)
	{
		if (value != (numOfProcesses - 1.0))
		{
			cout << "ERROR calculating ReduceToAverage double!!" << endl;
		}
		cout << "ReduceToAverage double done" << endl;
	}

	//ReduceToAverage double array test
	data[0] = processRank * 1.0; data[1] = processRank * 2.0; data[2] = processRank * 3.0;
	MPIMethods::ReduceToAverage(data, 3);
	if (processRank == rootRank)
	{
		if (data[0] != 0.5 * (numOfProcesses - 1.0) ||
			data[1] != 1.0 * (numOfProcesses - 1.0) ||
			data[2] != 1.5 * (numOfProcesses - 1.0))
		{
			cout << "ERROR calculating ReduceToAverage double array!!" << endl;
		}
		cout << "ReduceToAverage double array done" << endl;
	}

	//ReduceToAverage vector test
	v = {processRank * 1.0, processRank * 2.0, processRank * 3.0};
	MPIMethods::ReduceToAverage(v);
	if (processRank == rootRank)
	{
		if (v[0] != 0.5 * (numOfProcesses - 1.0) ||
			v[1] != 1.0 * (numOfProcesses - 1.0) ||
			v[2] != 1.5 * (numOfProcesses - 1.0))
		{
			cout << "ERROR calculating ReduceToAverage vector!!" << endl;
		}
		cout << "ReduceToAverage vector done" << endl;
	}

	//ReduceToAverage 2D vector test
	v2 = {{processRank * 1.0, processRank * 2.0, processRank * 3.0}, {processRank * 4.0, processRank * 5.0, processRank * 6.0}};
	MPIMethods::ReduceToAverage(v2);
	if (processRank == rootRank)
	{
		if (v2[0][0] != 0.5 * (numOfProcesses - 1.0) ||
			v2[0][1] != 1.0 * (numOfProcesses - 1.0) ||
			v2[0][2] != 1.5 * (numOfProcesses - 1.0) ||
			v2[1][0] != 2.0 * (numOfProcesses - 1.0) ||
			v2[1][1] != 2.5 * (numOfProcesses - 1.0) ||
			v2[1][2] != 3.0 * (numOfProcesses - 1.0))
		{
			cout << "ERROR calculating ReduceToAverage 2D vector!!" << endl;
		}
		cout << "ReduceToAverage 2D vector done" << endl;
	}

	//GatherForHistogram int test
	intValue = processRank % 3;
	h = MPIMethods::GatherHistogram(intValue);
	if (processRank == rootRank)
	{
		//TODO: check more properties of result
		if (h.size() != 3 ||
			h[0] != numOfProcesses / 3 + 1)
		{
			cout << "ERROR calculating GatherForHistogram int!!" << endl;
		}
		cout << "GatherForHistogram int done" << endl;
	}

	delete[] data;
}

void TestMPITiming(int multiplier, int nrOfRuns)
{
	Timer t;
	int k = nrOfRuns;
	int n = multiplier * multiplier;
	int n2D = multiplier;
	double* arr;
	vector<double> v;
	vector<vector<double> > v2D;

	if (processRank == rootRank)
	{
		cout << endl << "=====================" << endl << endl;
	}

	t.start();
	for (int i = 0; i < k; i++)
		CreateRandomVector(v, n, 0);
	t.stop();
	PrintTimings(MPIMethods::ReduceToMinMaxMean(t.duration()), "RandomVectorTest rand (" + to_string(n) + ")");

	t.start();
	for (int i = 0; i < k; i++)
		CreateRandomVector(v, n, 1);
	t.stop();
	PrintTimings(MPIMethods::ReduceToMinMaxMean(t.duration()), "RandomVectorTest mersenne (" + to_string(n) + ")");

	t.start();
	for (int i = 0; i < k; i++)
		CreateRandomVector(v, n, 2);
	t.stop();
	PrintTimings(MPIMethods::ReduceToMinMaxMean(t.duration()), "RandomVectorTest mersenne_64 (" + to_string(n) + ")");

	t.start();
	for (int i = 0; i < k; i++)
		CreateRandomVector2D(v2D, n2D, 0);
	t.stop();
	PrintTimings(MPIMethods::ReduceToMinMaxMean(t.duration()), "RandomVectorTest2D rand (" + to_string(n2D) + "x" + to_string(n2D) + ")");

	t.start();
	for (int i = 0; i < k; i++)
		CreateRandomVector2D(v2D, n2D, 1);
	t.stop();
	PrintTimings(MPIMethods::ReduceToMinMaxMean(t.duration()), "RandomVectorTest2D mersenne (" + to_string(n2D) + "x" + to_string(n2D) + ")");

	t.start();
	for (int i = 0; i < k; i++)
		CreateRandomVector2D(v2D, n2D, 2);
	t.stop();
	PrintTimings(MPIMethods::ReduceToMinMaxMean(t.duration()), "RandomVectorTest2D mersenne_64 (" + to_string(n2D) + "x" + to_string(n2D) + ")");

	if (processRank == rootRank)
	{
		cout << endl << "=====================" << endl << endl;
	}

	arr = new double[n];
	CreateRandomArray(arr, n);
	t.start();
	for (int i = 0; i < k; i++)
		MPIMethods::ReduceToAverage(arr, n);
	t.stop();
	PrintTimings(MPIMethods::ReduceToMinMaxMean(t.duration()), "MPIArrayTest (" + to_string(n) + ")");

	CreateRandomVector(v, n, 2);
	t.start();
	for (int i = 0; i < k; i++)
		MPIMethods::ReduceToAverage(v);
	t.stop();
	PrintTimings(MPIMethods::ReduceToMinMaxMean(t.duration()), "MPIVectorToArrayTest1 (" + to_string(n) + ")");

	//CreateRandomVector(v, n, 2);
	//t.start();
	//for (int i = 0; i < k; i++)
	//	MPIMethods::ReduceToAverage2(v);
	//t.stop();
	//PrintTimings(MPIMethods::ReduceToMinMaxMean(t.duration()), "MPIVectorToArrayTest2 (" + to_string(n) + ")");

	CreateRandomVector2D(v2D, n2D, 2);
	t.start();
	for (int i = 0; i < k; i++)
		MPIMethods::ReduceToAverage(v2D);
	t.stop();
	PrintTimings(MPIMethods::ReduceToMinMaxMean(t.duration()), "MPIVector2DToArrayTest (" + to_string(n2D) + "x" + to_string(n2D) + ")");

	delete[] arr;
}

void TestMPI(int argc, char *argv[])
{
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &processRank);
    MPI_Comm_size(MPI_COMM_WORLD, &numOfProcesses);
	//MPI::Init(argc, argv);
	//processRank = MPI::COMM_WORLD.Get_rank();
	//numOfProcesses = MPI::COMM_WORLD.Get_size();

	MPIMethods::numOfProcesses = numOfProcesses;
	MPIMethods::processRank = processRank;
	MPIMethods::rootRank = rootRank;


	if (processRank == rootRank)
	{
		cout << endl << "=====================";
		cout << endl << "=====================" << endl << endl;
	}
	GetCPUUsage();
	if (processRank == rootRank)
	{
		cout << endl << "==========TestMPIData()==========";
		cout << endl << "=================================" << endl << endl;
	}
	TestMPIData();
	if (processRank == rootRank)
	{
		cout << endl << "==========TestMPITiming(10, 100)==========";
		cout << endl << "==========================================" << endl << endl;
	}
	TestMPITiming(10, 100);
	if (processRank == rootRank)
	{
		cout << endl << "==========TestMPITiming(100, 100)==========";
		cout << endl << "===========================================" << endl << endl;
	}
	TestMPITiming(100, 100);
	if (processRank == rootRank)
	{
		cout << endl << "==========TestMPITiming(500, 100)==========";
		cout << endl << "===========================================" << endl << endl;
	}
	TestMPITiming(500, 100);
	if (processRank == rootRank)
	{
		cout << endl << "==========TestMPITiming(1000, 100)==========";
		cout << endl << "============================================" << endl << endl;
	}
	TestMPITiming(1000, 100);
	if (processRank == rootRank)
	{
		cout << endl << "=====================";
		cout << endl << "=====================" << endl << endl;
	}

	MPI_Finalize();
}

}

