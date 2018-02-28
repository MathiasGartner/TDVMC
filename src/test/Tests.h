#pragma once

#include "../Constants.h"
#include "../Utils.h"

using namespace std;

double eps = 1e-9;

bool CheckDouble(double d1, double d2)
{
	return fabs(d1 - d2) < eps;
}

void TestVectorDisplacement(double i1, double i2, double i3, double j1, double j2, double j3, double resultDisplacement, double result1, double result2, double result3)
{
	double* veci = new double[DIM];
	double* vecj = new double[DIM];
	double* vecij = new double[DIM];
	double vecDisplacement;

	veci[0] = i1;
	veci[1] = i2;
	veci[2] = i3;
	vecj[0] = j1;
	vecj[1] = j2;
	vecj[2] = j3;
	vecDisplacement = VectorDisplacementNIC(veci, vecj, vecij, DIM);
	if (CheckDouble(resultDisplacement, vecDisplacement) &&
		CheckDouble(result1, vecij[0]) &&
		CheckDouble(result2, vecij[1]) &&
		CheckDouble(result3, vecij[2]))
	{
		cout << "Test OK" << endl;
		cout << endl;
	}
	else
	{
		cout << "Test failed!!!!" << endl;
		cout << "vecDisplacement: " << resultDisplacement << "=" << vecDisplacement << endl;
		cout << "vecij[0]: " << result1 << "=" << vecij[0] << endl;
		cout << "vecij[0]: " << result2 << "=" << vecij[1] << endl;
		cout << "vecij[0]: " << result3 << "=" << vecij[2] << endl;
		cout << endl;
	}

	delete[] veci;
	delete[] vecj;
	delete[] vecij;
}

void TestVectorDisplacements()
{
	LBOX = 4;
	DIM = 3;
	TestVectorDisplacement(0, 0, 0, 0, 0, 0, 0, 0, 0, 0);
	TestVectorDisplacement(0, 0, 0, 1, 0, 0, 1, -1, 0, 0);
	TestVectorDisplacement(1, 0, 0, 0, 0, 0, 1, 1, 0, 0);
	TestVectorDisplacement(1, 0, 0, 6, 0, 0, 1, -1, 0, 0);
	TestVectorDisplacement(1, 0, 0, 3, 0, 0, 2, 2, 0, 0);
	TestVectorDisplacement(1, 0, 0, -1, 0, 0, 2, -2, 0, 0);
	TestVectorDisplacement(0, 0, 0, 10, 10, 10, sqrt(12), 2, 2, 2);
	TestVectorDisplacement(4, 4, 4, 10, 10, 10, sqrt(12), 2, 2, 2);
	TestVectorDisplacement(4, 4, 4, 10.1, 0, 0, sqrt(1.9*1.9), 1.9, 0, 0);
	TestVectorDisplacement(4, 4, 4, 9.9, 0, 0, sqrt(1.9*1.9), -1.9, 0, 0);
	TestVectorDisplacement(0.1, 0, 0, 3.9, 0, 0, sqrt(0.2*0.2), 0.2, 0, 0);
	cout << endl;
}
