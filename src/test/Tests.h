#pragma once

#include "../Constants.h"
#include "../Utils.h"

using namespace std;

namespace Test
{

double eps = 1e-9;

bool CheckDouble(double d1, double d2)
{
	return fabs(d1 - d2) < eps;
}

void TestVectorDisplacement(double i1, double i2, double i3, double j1, double j2, double j3, double resultDisplacement, double result1, double result2, double result3)
{
	vector<double> veci(DIM);
	vector<double> vecj(DIM);
	vector<double> vecij(DIM);
	double vecDisplacement;

	veci[0] = i1;
	veci[1] = i2;
	veci[2] = i3;
	vecj[0] = j1;
	vecj[1] = j2;
	vecj[2] = j3;
	vecDisplacement = VectorDisplacementNIC(veci, vecj, vecij);
	if (CheckDouble(resultDisplacement, vecDisplacement) &&
		CheckDouble(result1, vecij[0]) &&
		CheckDouble(result2, vecij[1]) &&
		CheckDouble(result3, vecij[2]))
	{
		cout << "Test OK" << endl;
	}
	else
	{
		cout << endl;
		cout << "Test failed!!!!" << endl;
		cout << "vecDisplacement: " << resultDisplacement << "=" << vecDisplacement << endl;
		cout << "vecij[0]: " << result1 << "=" << vecij[0] << endl;
		cout << "vecij[0]: " << result2 << "=" << vecij[1] << endl;
		cout << "vecij[0]: " << result3 << "=" << vecij[2] << endl;
		cout << endl;
	}
}

void TestVectorDisplacement(double i1, double j1, double resultDisplacement, double result1)
{
	vector<double> veci(DIM);
	vector<double> vecj(DIM);
	vector<double> vecij(DIM);
	double vecDisplacement;

	veci[0] = i1;
	vecj[0] = j1;
	vecDisplacement = VectorDisplacementNIC(veci, vecj, vecij);
	if (CheckDouble(resultDisplacement, vecDisplacement) &&
		CheckDouble(result1, vecij[0]))
	{
		cout << "Test OK" << endl;
	}
	else
	{
		cout << endl;
		cout << "Test failed!!!!" << endl;
		cout << "vecDisplacement: " << resultDisplacement << "=" << vecDisplacement << endl;
		cout << "vecij[0]: " << result1 << "=" << vecij[0] << endl;
		cout << endl;
	}
}

void TestVectorDisplacements()
{
	cout << "LBOX=4, DIM=3" << endl;
	LBOX = 4;
	LBOX_R = 1.0 / LBOX;
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
	TestVectorDisplacement(0, 0, 0, 3.9, 3.9, 3.9, sqrt(0.1*0.1*3), 0.1, 0.1, 0.1);

	TestVectorDisplacement(4.5, 2.1, -2.1, 0, 0, 0, sqrt(0.5*0.5 + 1.9*1.9 + 1.9*1.9), 0.5, -1.9, 1.9);
	TestVectorDisplacement(0.5, 1.5, 2.5, 0, 0, 0, sqrt(0.5*0.5 + 1.5*1.5 + 1.5*1.5), 0.5, 1.5, -1.5);
	TestVectorDisplacement(3.5, 4.5, 5.5, 0, 0, 0, sqrt(0.5*0.5 + 0.5*0.5 + 1.5*1.5), -0.5, 0.5, 1.5);
	cout << endl;

	cout << "LBOX=5, DIM=3" << endl;
	LBOX = 5;
	LBOX_R = 1.0 / LBOX;
	DIM = 3;
	TestVectorDisplacement(0, 0, 0, 0, 0, 0, 0, 0, 0, 0);
	TestVectorDisplacement(0, 0, 0, 1, 0, 0, 1, -1, 0, 0);
	TestVectorDisplacement(1, 0, 0, 0, 0, 0, 1, 1, 0, 0);
	TestVectorDisplacement(1, 0, 0, 6, 0, 0, 0, 0, 0, 0);
	TestVectorDisplacement(1, 0, 0, 3, 0, 0, 2, -2, 0, 0);
	TestVectorDisplacement(1, 0, 0, -1, 0, 0, 2, 2, 0, 0);
	TestVectorDisplacement(0, 0, 0, 10, 10, 10, 0, 0, 0, 0);
	TestVectorDisplacement(4, 4, 4, 10, 10, 10, sqrt(3), -1, -1, -1);
	TestVectorDisplacement(4, 4, 4, 10.1, 0, 0, sqrt(1.1*1.1 + 1 + 1), -1.1, -1, -1);
	TestVectorDisplacement(4, 4, 4, 9.9, 0, 0, sqrt(1 + 1 + 0.9*0.9), -0.9, -1, -1);
	TestVectorDisplacement(0.1, 0, 0, 3.9, 0, 0, 1.2, 1.2, 0, 0);
	TestVectorDisplacement(0, 0, 0, 4.9, 4.9, 4.9, sqrt(0.1*0.1*3), 0.1, 0.1, 0.1);

	TestVectorDisplacement(1.25, 0, 0, -1.25, 0, 0, 2.5, -2.5, 0, 0);
	TestVectorDisplacement(5.9, 2.6, -2.6, 0, 0, 0, sqrt(0.9*0.9 + 2.4*2.4 + 2.4*2.4), 0.9, -2.4, 2.4);

	TestVectorDisplacement(0, 1, 2, 0, 0, 0, sqrt(1 + 2*2), 0, 1, 2);
	TestVectorDisplacement(3, 4, 5, 0, 0, 0, sqrt(2*2 + 1), -2, -1, 0);
	cout << endl;

	cout << "LBOX=4, DIM=1" << endl;
	LBOX = 4;
	LBOX_R = 1.0 / LBOX;
	DIM = 1;
	TestVectorDisplacement(0, 0, 0, 0);
	TestVectorDisplacement(0, 1, 1, -1);
	TestVectorDisplacement(1, 0, 1, 1);
	TestVectorDisplacement(1, 1, 0, 0);
	TestVectorDisplacement(0, 2, 2, 2); //INFO: or -2?
	TestVectorDisplacement(0, 3, 1, 1);
	TestVectorDisplacement(0, 4, 0, 0);
	TestVectorDisplacement(0, 5, 1, -1);
	cout << endl;

	cout << "LBOX=5, DIM=1" << endl;
	LBOX = 5;
	LBOX_R = 1.0 / LBOX;
	DIM = 1;
	TestVectorDisplacement(0, 0, 0, 0);
	TestVectorDisplacement(0, 1, 1, -1);
	TestVectorDisplacement(1, 0, 1, 1);
	TestVectorDisplacement(1, 1, 0, 0);
	TestVectorDisplacement(0, 2, 2, -2);
	TestVectorDisplacement(0, 3, 2, 2);
	TestVectorDisplacement(0, 4, 1, 1);
	TestVectorDisplacement(0, 5, 0, 0);
	TestVectorDisplacement(0, 6, 1, -1);
	cout << endl;
}

}
