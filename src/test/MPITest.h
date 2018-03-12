#pragma once

#include "../Constants.h"
#include "../Timer.h"
#include "../Utils.h"
#include "mpi.h"

using namespace std;

int ownRank;
int totalSize;

void MPIVectorToArrayTest(vector<double>& v)
{
	int rootRank = 0;

	double* ownValues = new double[v.size()];
	double* reducedValues = new double[v.size()];

	for (unsigned int i = 0; i < v.size(); i++)
	{
		ownValues[i] = v[i];
	}
	MPI_Reduce(ownValues, reducedValues, v.size(), MPI_DOUBLE, MPI_SUM, rootRank, MPI_COMM_WORLD);
	if (ownRank == rootRank)
	{
		for (unsigned int i = 0; i < v.size(); i++)
		{
			v[i] = reducedValues[i] / (double)totalSize;
		}
	}
	delete[] ownValues;
	delete[] reducedValues;
}

void MPIVector2DToArrayTest(vector<vector<double> >& v)
{
	for (unsigned int i = 0; i < v.size(); i++)
	{
		MPIVectorToArrayTest(v[i]);
	}
}

vector<double> GetTimings(double* data)
{
	int rootRank = 0;
	vector<double> values = {};
	double own = 0;
	double min = 0;
	double max = 0;
	double sum = 0;

	own = *data;
	MPI_Reduce(&own, &min, 1, MPI_DOUBLE, MPI_MIN, rootRank, MPI_COMM_WORLD);
	MPI_Reduce(&own, &max, 1, MPI_DOUBLE, MPI_MAX, rootRank, MPI_COMM_WORLD);
	MPI_Reduce(&own, &sum, 1, MPI_DOUBLE, MPI_SUM, rootRank, MPI_COMM_WORLD);
	if (ownRank == rootRank)
	{
		values = { min, max, sum / (double)totalSize };
	}
	return values;
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
	}
	else if (randType == 2)
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
	int rootRank = 0;
	if (ownRank == rootRank)
	{
		cout << message << endl;
		cout << "duration: min = " << to_string(timings[0]) << " ms" << endl;
		cout << "          max = " << to_string(timings[1]) << " ms" << endl;
		cout << "          <t> = " << to_string(timings[2]) << " ms" << endl;
	}
}

void TestMPI(int argc, char *argv[])
{
	MPI::Init(argc, argv);
	ownRank = MPI::COMM_WORLD.Get_rank();
	totalSize = MPI::COMM_WORLD.Get_size();

	Timer t;
	double duration;
	vector<double> timings;
	int k = 1;
	int n = 1000000;
	int n2D = 1000;
	vector<double> v;
	vector<vector<double> > v2D;

	t.start();
	for (int i = 0; i < k; i++)
	CreateRandomVector(v, n, 0);
	t.stop();
	//cout << "RandomVectorTest rand with n=" << to_string(n) << ": " << to_string(t.duration()) << "ms" << endl;
	duration = t.duration();
	timings = GetTimings(&duration);
	PrintTimings(timings, "RandomVectorTest rand");

	t.start();
	for (int i = 0; i < k; i++)
	CreateRandomVector(v, n, 1);
	t.stop();
	//cout << "RandomVectorTest mersenne with n=" << to_string(n) << ": " << to_string(t.duration()) << "ms" << endl;
	duration = t.duration();
	timings = GetTimings(&duration);
	PrintTimings(timings, "RandomVectorTest mersenne");

	t.start();
	for (int i = 0; i < k; i++)
	CreateRandomVector(v, n, 2);
	t.stop();
	//cout << "RandomVectorTest mersenne_64 with n=" << to_string(n) << ": " << to_string(t.duration()) << "ms" << endl;
	duration = t.duration();
	timings = GetTimings(&duration);
	PrintTimings(timings, "RandomVectorTest mersenne_64");

	t.start();
	for (int i = 0; i < k; i++)
	CreateRandomVector2D(v2D, n2D, 0);
	t.stop();
	//cout << "RandomVectorTest2D rand with n=" << to_string(n) << ": " << to_string(t.duration()) << "ms" << endl;
	duration = t.duration();
	timings = GetTimings(&duration);
	PrintTimings(timings, "RandomVectorTest2D rand");

	t.start();
	for (int i = 0; i < k; i++)
	CreateRandomVector2D(v2D, n2D, 1);
	t.stop();
	//cout << "RandomVectorTest2D mersenne with n=" << to_string(n) << ": " << to_string(t.duration()) << "ms" << endl;
	duration = t.duration();
	timings = GetTimings(&duration);
	PrintTimings(timings, "RandomVectorTest2D rand");

	t.start();
	for (int i = 0; i < k; i++)
	CreateRandomVector2D(v2D, n2D, 2);
	t.stop();
	//cout << "RandomVectorTest2D mersenne_64 with n=" << to_string(n) << ": " << to_string(t.duration()) << "ms" << endl;
	duration = t.duration();
	timings = GetTimings(&duration);
	PrintTimings(timings, "RandomVectorTest2D rand");

	t.start();
	for (int i = 0; i < k; i++)
	MPIVectorToArrayTest(v);
	t.stop();
	//cout << "MPIVectorToArrayTest with n=" << to_string(n) << ": " << to_string(t.duration()) << "ms" << endl;
	duration = t.duration();
	timings = GetTimings(&duration);
	PrintTimings(timings, "MPIVectorToArrayTest");

	t.start();
	for (int i = 0; i < k; i++)
	MPIVector2DToArrayTest(v2D);
	t.stop();
	//cout << "MPIVector2DToArrayTest with n=" << to_string(n2D) << ": " << to_string(t.duration()) << "ms" << endl;
	duration = t.duration();
	timings = GetTimings(&duration);
	PrintTimings(timings, "MPIVector2DToArrayTest");

	MPI::Finalize();
}
