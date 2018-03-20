#pragma once

#include <math.h>
#include "mpi.h"
#include <iostream>
#include <stdlib.h>

using namespace std;

namespace Test
{

void CalculatePi(int argc, char *argv[]) {
	int rank;
	int size;
	double PI25DT = 3.141592653589793238462643;
	double mypi, pi, x, y;
	int total, inside;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

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

	MPI_Reduce(&mypi, &pi, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	pi = pi / (double)size;
	if (rank == 0)
		cout << "pi is approximately " << pi << ", Error is "
				<< fabs(pi - PI25DT) << endl;

	MPI_Finalize();
}

}
