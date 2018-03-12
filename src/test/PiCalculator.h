#pragma once

#include <math.h>
#include "mpi.h"
#include <iostream>
#include <stdlib.h>

using namespace std;

void CalculatePi(int argc, char *argv[]) {
	int rank;
	int size;
	double PI25DT = 3.141592653589793238462643;
	double mypi, pi, x, y;
	int total, inside;

	MPI::Init(argc, argv);
	size = MPI::COMM_WORLD.Get_size();
	rank = MPI::COMM_WORLD.Get_rank();

	total = 10000;
	srand(rank + 1);

	inside = 0;
	for (int i = 0; i < total; i++)
	{
		x = ((double) rand() / RAND_MAX);
		y = ((double) rand() / RAND_MAX);
		if (sqrt(x*x + y*y) < 1.0)
		{
			inside++;
		}
	}
	mypi = ((double)inside) / ((double)total) * 4.0;
	cout << "this is process# " << rank << " mypi=" << mypi << endl;

	MPI::COMM_WORLD.Reduce(&mypi, &pi, 1, MPI::DOUBLE, MPI::SUM, 0);
	pi = pi / (double)size;
	if (rank == 0)
		cout << "pi is approximately " << pi << ", Error is "
				<< fabs(pi - PI25DT) << endl;

	MPI::Finalize();
}
