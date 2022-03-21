#include "SplineFactory.h"

using namespace std;

namespace SplineFactory
{

vector<vector<vector<double> > > GetWeights(vector<double> n)
{
	return GetWeights3(n);
}

vector<vector<vector<double> > > GetWeights(vector<double> n, int splineOrder)
{
	if (splineOrder == 0)
	{
		return GetWeights0(n);
	}
	else if (splineOrder == 3)
	{
		return GetWeights3(n);
	}
	else if (splineOrder == 4)
	{
		return GetWeights4(n);
	}
	else
	{
		return GetWeights3(n);
	}
}

vector<vector<vector<double> > > GetWeights0(vector<double> n)
{
	vector<vector<vector<double> > > weights; //INFO: indices (s, p, c): s...spline, p...part of spline, c...polynomial coefficients

	int nOrder = 0;
	int nParts;
	int nNodes;
	int nSplines;

	nParts = nOrder + 1;
	nNodes = n.size();
	nSplines = nNodes - nParts;
	InitVector(weights, nSplines, nParts, nOrder + 1, 0.0);

	for (int s = 0; s < nSplines; s++)
	{
		weights[s][0][0] = 1.0;
	}
	return weights;
}

vector<vector<vector<double> > > GetWeights3(vector<double> n)
{
	//see file NonUniformSplines_new_v2.nb -> Own spline implementation -> third order -> Part 1/2/3/4
	vector<vector<vector<double> > > weights; //INFO: indices (s, p, c): s...spline, p...part of spline, c...polynomial coefficients

	int nOrder = 3;
	int nParts;
	int nNodes;
	int nSplines;

	nParts = nOrder + 1;
	nNodes = n.size();
	nSplines = nNodes - nParts;
	weights.resize(nSplines);
	for (auto& w : weights)
	{
		w.resize(nParts);
		for (auto& p : w)
		{
			p.resize(nOrder + 1);
		}
	}

	for (int s = 0; s < nSplines; s++)
	{
		int i = s;
		//Part 1:
		weights[s][0][0] = pow(n[i], 3) / ((n[i] - n[i + 1]) * (n[i] - n[i + 2]) * (n[i] - n[i + 3]));
		weights[s][0][1] = (-3 * pow(n[i], 2)) / ((n[i] - n[i + 1]) * (n[i] - n[i + 2]) * (n[i] - n[i + 3]));
		weights[s][0][2] = (3 * n[i]) / ((n[i] - n[i + 1]) * (n[i] - n[i + 2]) * (n[i] - n[i + 3]));
		weights[s][0][3] = -(1 / ((n[i] - n[i + 1]) * (n[i] - n[i + 2]) * (n[i] - n[i + 3])));
		//Part 2:
		weights[s][1][0] = (pow(n[i], 2) * pow(n[i + 1], 2) * n[i + 2] + pow(n[i], 2) * pow(n[i + 1], 2) * n[i + 3] + pow(n[i], 2) * pow(n[i + 1], 2) * n[i + 4] - pow(n[i], 2) * n[i + 1] * n[i + 3] * n[i + 4] - n[i] * pow(n[i + 1], 2) * n[i + 3] * n[i + 4] + pow(n[i], 2) * n[i + 2] * n[i + 3] * n[i + 4] + n[i] * n[i + 1] * n[i + 2] * n[i + 3] * n[i + 4] + pow(n[i + 1], 2) * n[i + 2] * n[i + 3] * n[i + 4] - pow(n[i], 2) * n[i + 1] * n[i + 2] * (n[i + 3] + n[i + 4]) - n[i] * pow(n[i + 1], 2) * n[i + 2] * (n[i + 3] + n[i + 4])) / ((n[i + 1] - n[i + 2]) * (-n[i] + n[i + 2]) * (n[i] - n[i + 3]) * (n[i + 1] - n[i + 3]) * (n[i + 1] - n[i + 4]));
		weights[s][1][1] = (-3 * pow(n[i], 2) * pow(n[i + 1], 2) + 3 * n[i] * n[i + 1] * n[i + 2] * n[i + 3] + 3 * n[i] * n[i + 1] * n[i + 2] * n[i + 4] + 3 * n[i] * n[i + 1] * n[i + 3] * n[i + 4] - 3 * n[i] * n[i + 2] * n[i + 3] * n[i + 4] - 3 * n[i + 1] * n[i + 2] * n[i + 3] * n[i + 4]) / ((n[i + 1] - n[i + 2]) * (-n[i] + n[i + 2]) * (n[i] - n[i + 3]) * (n[i + 1] - n[i + 3]) * (n[i + 1] - n[i + 4]));
		weights[s][1][2] = (3 * pow(n[i], 2) * n[i + 1] + 3 * n[i] * pow(n[i + 1], 2) - 3 * n[i] * n[i + 1] * n[i + 2] - 3 * n[i] * n[i + 1] * n[i + 3] - 3 * n[i] * n[i + 1] * n[i + 4] + 3 * n[i + 2] * n[i + 3] * n[i + 4]) / ((n[i + 1] - n[i + 2]) * (-n[i] + n[i + 2]) * (n[i] - n[i + 3]) * (n[i + 1] - n[i + 3]) * (n[i + 1] - n[i + 4]));
		weights[s][1][3] = (-pow(n[i], 2) - n[i] * n[i + 1] - pow(n[i + 1], 2) + n[i] * n[i + 2] + n[i + 1] * n[i + 2] - n[i + 2] * n[i + 3] - n[i + 2] * n[i + 4] - n[i + 3] * n[i + 4] + n[i] * (n[i + 3] + n[i + 4]) + n[i + 1] * (n[i + 3] + n[i + 4])) / ((n[i + 1] - n[i + 2]) * (-n[i] + n[i + 2]) * (n[i] - n[i + 3]) * (n[i + 1] - n[i + 3]) * (n[i + 1] - n[i + 4]));
		//Part 3:
		weights[s][2][0] = (n[i] * n[i + 1] * n[i + 2] * pow(n[i + 3], 2) + n[i] * n[i + 1] * n[i + 2] * n[i + 3] * n[i + 4] - n[i] * n[i + 1] * pow(n[i + 3], 2) * n[i + 4] - n[i] * n[i + 2] * pow(n[i + 3], 2) * n[i + 4] - n[i + 1] * n[i + 2] * pow(n[i + 3], 2) * n[i + 4] + n[i] * n[i + 1] * n[i + 2] * pow(n[i + 4], 2) - n[i] * n[i + 1] * n[i + 3] * pow(n[i + 4], 2) - n[i] * n[i + 2] * n[i + 3] * pow(n[i + 4], 2) - n[i + 1] * n[i + 2] * n[i + 3] * pow(n[i + 4], 2) + n[i] * pow(n[i + 3], 2) * pow(n[i + 4], 2) + n[i + 1] * pow(n[i + 3], 2) * pow(n[i + 4], 2) + n[i + 2] * pow(n[i + 3], 2) * pow(n[i + 4], 2)) / ((n[i + 2] - n[i + 3]) * (-n[i] + n[i + 3]) * (-n[i + 1] + n[i + 3]) * (n[i + 1] - n[i + 4]) * (n[i + 2] - n[i + 4]));
		weights[s][2][1] = (-3 * n[i] * n[i + 1] * n[i + 2] * n[i + 3] - 3 * n[i] * n[i + 1] * n[i + 2] * n[i + 4] + 3 * n[i] * n[i + 1] * n[i + 3] * n[i + 4] + 3 * n[i] * n[i + 2] * n[i + 3] * n[i + 4] + 3 * n[i + 1] * n[i + 2] * n[i + 3] * n[i + 4] - 3 * pow(n[i + 3], 2) * pow(n[i + 4], 2)) / ((n[i + 2] - n[i + 3]) * (-n[i] + n[i + 3]) * (-n[i + 1] + n[i + 3]) * (n[i + 1] - n[i + 4]) * (n[i + 2] - n[i + 4]));
		weights[s][2][2] = (3 * n[i] * n[i + 1] * n[i + 2] - 3 * n[i] * n[i + 3] * n[i + 4] - 3 * n[i + 1] * n[i + 3] * n[i + 4] - 3 * n[i + 2] * n[i + 3] * n[i + 4] + 3 * pow(n[i + 3], 2) * n[i + 4] + 3 * n[i + 3] * pow(n[i + 4], 2)) / ((n[i + 2] - n[i + 3]) * (-n[i] + n[i + 3]) * (-n[i + 1] + n[i + 3]) * (n[i + 1] - n[i + 4]) * (n[i + 2] - n[i + 4]));
		weights[s][2][3] = (-(n[i] * n[i + 1]) - n[i] * n[i + 2] - n[i + 1] * n[i + 2] + n[i] * n[i + 3] + n[i + 1] * n[i + 3] + n[i + 2] * n[i + 3] - pow(n[i + 3], 2) + n[i] * n[i + 4] + n[i + 1] * n[i + 4] + n[i + 2] * n[i + 4] - n[i + 3] * n[i + 4] - pow(n[i + 4], 2)) / ((n[i + 2] - n[i + 3]) * (-n[i] + n[i + 3]) * (-n[i + 1] + n[i + 3]) * (n[i + 1] - n[i + 4]) * (n[i + 2] - n[i + 4]));
		//Part 4:
		weights[s][3][0] = -(pow(n[i + 4], 3) / ((n[i + 3] - n[i + 4]) * (-n[i + 1] + n[i + 4]) * (-n[i + 2] + n[i + 4])));
		weights[s][3][1] = (3 * pow(n[i + 4], 2)) / ((n[i + 3] - n[i + 4]) * (-n[i + 1] + n[i + 4]) * (-n[i + 2] + n[i + 4]));
		weights[s][3][2] = (-3 * n[i + 4]) / ((n[i + 3] - n[i + 4]) * (-n[i + 1] + n[i + 4]) * (-n[i + 2] + n[i + 4]));
		weights[s][3][3] = 1 / ((n[i + 3] - n[i + 4]) * (-n[i + 1] + n[i + 4]) * (-n[i + 2] + n[i + 4]));
	}
	return weights;
}

vector<vector<vector<double> > > GetWeights4(vector<double> n)
{
	vector<vector<vector<double> > > weights; //INFO: indices (s, p, c): s...spline, p...part of spline, c...polynomial coefficients

	int nOrder = 4;
	int nParts;
	int nNodes;
	int nSplines;

	nParts = nOrder + 1;
	nNodes = n.size();
	nSplines = nNodes - nParts;
	InitVector(weights, nSplines, nParts, nOrder + 1, 0.0);

	for (int s = 0; s < nSplines; s++)
	{
		int i = s;
		//Part 1:
		weights[s][0][0] = pow(n[i], 4) / ((n[i] - n[i + 1]) * (n[i] - n[i + 2]) * (n[i] - n[i + 3]) * (n[i] - n[i + 4]));
		weights[s][0][1] = (-4 * pow(n[i], 3)) / ((n[i] - n[i + 1]) * (n[i] - n[i + 2]) * (n[i] - n[i + 3]) * (n[i] - n[i + 4]));
		weights[s][0][2] = (6 * pow(n[i], 2)) / ((n[i] - n[i + 1]) * (n[i] - n[i + 2]) * (n[i] - n[i + 3]) * (n[i] - n[i + 4]));
		weights[s][0][3] = (-4.0 * n[i]) / ((n[i] - n[i + 1]) * (n[i] - n[i + 2]) * (n[i] - n[i + 3]) * (n[i] - n[i + 4]));
		weights[s][0][4] = 1.0 / ((n[i] - n[i + 1]) * (n[i] - n[i + 2]) * (n[i] - n[i + 3]) * (n[i] - n[i + 4]));
		//Part 2:
		weights[s][1][0] = (-(pow(n[i + 1], 3) * n[i + 2] * n[i + 3] * n[i + 4] * n[i + 5]) + pow(n[i], 3) * (n[i + 1] * (n[i + 2] * n[i + 3] * n[i + 4] + pow(n[i + 1], 2) * (n[i + 2] + n[i + 3] + n[i + 4]) - n[i + 1] * (n[i + 2] * n[i + 3] + (n[i + 2] + n[i + 3]) * n[i + 4])) + (n[i + 1] - n[i + 2]) * (n[i + 1] - n[i + 3]) * (n[i + 1] - n[i + 4]) * n[i + 5]) + n[i] * pow(n[i + 1], 2) * (n[i + 1] * n[i + 2] * n[i + 3] * n[i + 4] + (n[i + 1] * n[i + 2] * n[i + 3] - n[i + 2] * n[i + 3] * n[i + 4] + n[i + 1] * (n[i + 2] + n[i + 3]) * n[i + 4]) * n[i + 5])
				+ pow(n[i], 2) * n[i + 1] * (-(n[i + 2] * n[i + 3] * n[i + 4] * n[i + 5]) + n[i + 1] * (n[i + 2] * n[i + 3] * n[i + 4] + n[i + 3] * n[i + 4] * n[i + 5] + n[i + 2] * (n[i + 3] + n[i + 4]) * n[i + 5]) - pow(n[i + 1], 2) * (n[i + 3] * n[i + 4] + (n[i + 3] + n[i + 4]) * n[i + 5] + n[i + 2] * (n[i + 3] + n[i + 4] + n[i + 5])))) / ((n[i + 1] - n[i + 2]) * (-n[i] + n[i + 2]) * (n[i] - n[i + 3]) * (n[i + 1] - n[i + 3]) * (n[i] - n[i + 4]) * (n[i + 1] - n[i + 4]) * (n[i + 1] - n[i + 5]));
		weights[s][1][1] = (-4 * pow(n[i], 3) * pow(n[i + 1], 3) + 4 * pow(n[i + 1], 2) * n[i + 2] * n[i + 3] * n[i + 4] * n[i + 5] - 4 * n[i] * n[i + 1] * (n[i + 1] * n[i + 2] * n[i + 3] * n[i + 4] + (n[i + 1] * n[i + 2] * n[i + 3] - n[i + 2] * n[i + 3] * n[i + 4] + n[i + 1] * (n[i + 2] + n[i + 3]) * n[i + 4]) * n[i + 5]) + 4 * pow(n[i], 2) * (n[i + 2] * n[i + 3] * n[i + 4] * n[i + 5] - n[i + 1] * (n[i + 2] * n[i + 3] * n[i + 4] + n[i + 3] * n[i + 4] * n[i + 5] + n[i + 2] * (n[i + 3] + n[i + 4]) * n[i + 5]) + pow(n[i + 1], 2) * (n[i + 3] * n[i + 4] + (n[i + 3] + n[i + 4]) * n[i + 5] + n[i + 2] * (n[i + 3] + n[i + 4] + n[i + 5])))) / ((n[i + 1] - n[i + 2]) * (-n[i] + n[i + 2]) * (n[i] - n[i + 3]) * (n[i + 1] - n[i + 3]) * (n[i] - n[i + 4]) * (n[i + 1] - n[i + 4]) * (n[i + 1] - n[i + 5]));
		weights[s][1][2] = (6 * n[i] * n[i + 1] * (n[i] * n[i + 1] * (n[i] + n[i + 1] - n[i + 2] - n[i + 3]) - n[i] * n[i + 1] * n[i + 4] + n[i + 2] * n[i + 3] * n[i + 4]) + 6 * (n[i] * n[i + 1] * (-(n[i] * n[i + 1]) + n[i + 2] * n[i + 3]) + (-(n[i + 1] * n[i + 2] * n[i + 3]) + n[i] * (-(n[i + 2] * n[i + 3]) + n[i + 1] * (n[i + 2] + n[i + 3]))) * n[i + 4]) * n[i + 5]) / ((n[i + 1] - n[i + 2]) * (-n[i] + n[i + 2]) * (n[i] - n[i + 3]) * (n[i + 1] - n[i + 3]) * (n[i] - n[i + 4]) * (n[i + 1] - n[i + 4]) * (n[i + 1] - n[i + 5]));
		weights[s][1][3] = (-4 * pow(n[i], 3) * n[i + 1] + 4 * n[i + 2] * n[i + 3] * n[i + 4] * n[i + 5] + 4 * pow(n[i], 2) * n[i + 1] * (-n[i + 1] + n[i + 2] + n[i + 3] + n[i + 4] + n[i + 5]) - 4 * n[i] * n[i + 1] * (pow(n[i + 1], 2) + n[i + 3] * n[i + 4] + (n[i + 3] + n[i + 4]) * n[i + 5] + n[i + 2] * (n[i + 3] + n[i + 4] + n[i + 5]) - n[i + 1] * (n[i + 2] + n[i + 3] + n[i + 4] + n[i + 5]))) / ((n[i + 1] - n[i + 2]) * (-n[i] + n[i + 2]) * (n[i] - n[i + 3]) * (n[i + 1] - n[i + 3]) * (n[i] - n[i + 4]) * (n[i + 1] - n[i + 4]) * (n[i + 1] - n[i + 5]));
		weights[s][1][4] = (-((pow(n[i], 2) + pow(n[i + 1], 2) + n[i + 2] * n[i + 3] + (n[i + 2] + n[i + 3]) * n[i + 4] - n[i + 1] * (n[i + 2] + n[i + 3] + n[i + 4]) - n[i] * (-n[i + 1] + n[i + 2] + n[i + 3] + n[i + 4])) / ((n[i] - n[i + 2]) * (n[i] - n[i + 3]) * (n[i] - n[i + 4]))) + 1 / (-n[i + 1] + n[i + 5])) / ((n[i + 1] - n[i + 2]) * (n[i + 1] - n[i + 3]) * (n[i + 1] - n[i + 4]));
		//Part 3:
		weights[s][2][0] = (n[i + 3] * n[i + 4] * n[i + 5] * (pow(n[i + 2], 2) * n[i + 3] * n[i + 4] * n[i + 5] - n[i + 1] * n[i + 2] * (n[i + 2] * n[i + 3] * n[i + 4] - n[i + 3] * n[i + 4] * n[i + 5] + n[i + 2] * (n[i + 3] + n[i + 4]) * n[i + 5]) + pow(n[i + 1], 2) * (n[i + 3] * n[i + 4] * n[i + 5] + pow(n[i + 2], 2) * (n[i + 3] + n[i + 4] + n[i + 5]) - n[i + 2] * (n[i + 3] * n[i + 4] + (n[i + 3] + n[i + 4]) * n[i + 5])))
				- n[i] * (n[i + 2] * n[i + 3] * n[i + 4] * n[i + 5] * (n[i + 2] * n[i + 3] * n[i + 4] - n[i + 3] * n[i + 4] * n[i + 5] + n[i + 2] * (n[i + 3] + n[i + 4]) * n[i + 5]) - n[i + 1] * pow(n[i + 2] * n[i + 3] * n[i + 4] - n[i + 3] * n[i + 4] * n[i + 5] + n[i + 2] * (n[i + 3] + n[i + 4]) * n[i + 5], 2) + pow(n[i + 1], 2) * (pow(n[i + 2], 2) * (n[i + 3] + n[i + 4]) * (n[i + 3] + n[i + 5]) * (n[i + 4] + n[i + 5]) + n[i + 3] * n[i + 4] * n[i + 5] * (n[i + 3] * n[i + 4] + (n[i + 3] + n[i + 4]) * n[i + 5]) - n[i + 2] * pow(n[i + 3] * n[i + 4] + (n[i + 3] + n[i + 4]) * n[i + 5], 2)))
				+ pow(n[i], 2)
						* (n[i + 3] * n[i + 4] * n[i + 5] * (n[i + 3] * n[i + 4] * n[i + 5] + pow(n[i + 2], 2) * (n[i + 3] + n[i + 4] + n[i + 5]) - n[i + 2] * (n[i + 3] * n[i + 4] + (n[i + 3] + n[i + 4]) * n[i + 5])) + n[i + 1] * (-(pow(n[i + 2], 2) * (n[i + 3] + n[i + 4]) * (n[i + 3] + n[i + 5]) * (n[i + 4] + n[i + 5])) - n[i + 3] * n[i + 4] * n[i + 5] * (n[i + 3] * n[i + 4] + (n[i + 3] + n[i + 4]) * n[i + 5]) + n[i + 2] * pow(n[i + 3] * n[i + 4] + (n[i + 3] + n[i + 4]) * n[i + 5], 2))
								+ pow(n[i + 1], 2) * (-(n[i + 2] * (n[i + 3] + n[i + 4]) * (n[i + 3] + n[i + 5]) * (n[i + 4] + n[i + 5])) + n[i + 3] * n[i + 4] * n[i + 5] * (n[i + 3] + n[i + 4] + n[i + 5]) + pow(n[i + 2], 2) * (pow(n[i + 3], 2) + n[i + 3] * n[i + 4] + pow(n[i + 4], 2) + (n[i + 3] + n[i + 4]) * n[i + 5] + pow(n[i + 5], 2))))) / ((n[i + 2] - n[i + 3]) * (-n[i] + n[i + 3]) * (-n[i + 1] + n[i + 3]) * (n[i] - n[i + 4]) * (n[i + 1] - n[i + 4]) * (n[i + 2] - n[i + 4]) * (n[i + 1] - n[i + 5]) * (n[i + 2] - n[i + 5]));
		weights[s][2][1] = (-4 * n[i] * n[i + 1] * n[i + 2] * pow(n[i + 3], 2) + 4 * n[i + 3] * (-(n[i] * n[i + 1] * n[i + 2]) + n[i + 1] * n[i + 2] * n[i + 3] + n[i] * (n[i + 1] + n[i + 2]) * n[i + 3]) * n[i + 4] - 4 * (n[i] * n[i + 1] * n[i + 2] - (n[i + 1] * n[i + 2] + n[i] * (n[i + 1] + n[i + 2])) * n[i + 3] + (n[i] + n[i + 1] + n[i + 2]) * pow(n[i + 3], 2)) * pow(n[i + 4], 2)) / ((n[i + 2] - n[i + 3]) * (-n[i] + n[i + 3]) * (-n[i + 1] + n[i + 3]) * (n[i] - n[i + 4]) * (n[i + 1] - n[i + 4]) * (n[i + 2] - n[i + 4])) + (4 * pow(n[i + 1], 3)) / ((n[i + 1] - n[i + 2]) * (n[i + 1] - n[i + 3]) * (n[i + 1] - n[i + 4]) * (n[i + 1] - n[i + 5])) + (4 * pow(n[i + 2], 3)) / ((-n[i + 1] + n[i + 2]) * (n[i + 2] - n[i + 3]) * (n[i + 2] - n[i + 4]) * (n[i + 2] - n[i + 5]));
		weights[s][2][2] = 6 * ((n[i] * n[i + 1] * n[i + 2] * n[i + 3] + (n[i] * n[i + 1] * n[i + 2] - (n[i + 1] * n[i + 2] + n[i] * (n[i + 1] + n[i + 2])) * n[i + 3]) * n[i + 4] + pow(n[i + 3], 2) * pow(n[i + 4], 2)) / ((n[i + 2] - n[i + 3]) * (-n[i] + n[i + 3]) * (-n[i + 1] + n[i + 3]) * (n[i] - n[i + 4]) * (n[i + 1] - n[i + 4]) * (n[i + 2] - n[i + 4])) + pow(n[i + 2], 2) / ((n[i + 1] - n[i + 2]) * (n[i + 2] - n[i + 3]) * (n[i + 2] - n[i + 4]) * (n[i + 2] - n[i + 5])) + pow(n[i + 1], 2) / ((n[i + 1] - n[i + 2]) * (n[i + 1] - n[i + 3]) * (n[i + 1] - n[i + 4]) * (-n[i + 1] + n[i + 5])));
		weights[s][2][3] = (-4 * n[i] * n[i + 1] * n[i + 2] + 4 * (n[i] + n[i + 1] + n[i + 2] - n[i + 3]) * n[i + 3] * n[i + 4] - 4 * n[i + 3] * pow(n[i + 4], 2)) / ((n[i + 2] - n[i + 3]) * (-n[i] + n[i + 3]) * (-n[i + 1] + n[i + 3]) * (n[i] - n[i + 4]) * (n[i + 1] - n[i + 4]) * (n[i + 2] - n[i + 4])) + (4 * n[i + 1]) / ((n[i + 1] - n[i + 2]) * (n[i + 1] - n[i + 3]) * (n[i + 1] - n[i + 4]) * (n[i + 1] - n[i + 5])) + (4 * n[i + 2]) / ((-n[i + 1] + n[i + 2]) * (n[i + 2] - n[i + 3]) * (n[i + 2] - n[i + 4]) * (n[i + 2] - n[i + 5]));
		weights[s][2][4] = ((n[i + 1] - n[i + 3]) * (n[i + 2] - n[i + 3]) + n[i] * (n[i + 1] + n[i + 2] - n[i + 3] - n[i + 4]) - (n[i + 1] + n[i + 2] - n[i + 3]) * n[i + 4] + pow(n[i + 4], 2)) / ((n[i + 2] - n[i + 3]) * (-n[i] + n[i + 3]) * (-n[i + 1] + n[i + 3]) * (n[i] - n[i + 4]) * (n[i + 1] - n[i + 4]) * (n[i + 2] - n[i + 4])) + 1 / ((-n[i + 1] + n[i + 2]) * (n[i + 1] - n[i + 3]) * (n[i + 1] - n[i + 4]) * (n[i + 1] - n[i + 5])) + 1 / ((n[i + 1] - n[i + 2]) * (n[i + 2] - n[i + 3]) * (n[i + 2] - n[i + 4]) * (n[i + 2] - n[i + 5]));
		//Part 4:
		weights[s][3][0] = (n[i] * n[i + 1] * n[i + 2] * n[i + 3] * pow(n[i + 4], 3) - pow(n[i + 4], 2) * (n[i + 1] * n[i + 2] * n[i + 3] * n[i + 4] + n[i] * (-(n[i + 1] * n[i + 2] * n[i + 3]) + n[i + 2] * n[i + 3] * n[i + 4] + n[i + 1] * (n[i + 2] + n[i + 3]) * n[i + 4])) * n[i + 5] + n[i + 4] * (n[i] * n[i + 1] * n[i + 2] * n[i + 3] - (n[i] * n[i + 1] * n[i + 2] + n[i + 1] * n[i + 2] * n[i + 3] + n[i] * (n[i + 1] + n[i + 2]) * n[i + 3]) * n[i + 4] + (n[i + 2] * n[i + 3] + n[i + 1] * (n[i + 2] + n[i + 3]) + n[i] * (n[i + 1] + n[i + 2] + n[i + 3])) * pow(n[i + 4], 2)) * pow(n[i + 5], 2)
				+ (n[i] * (n[i + 1] - n[i + 4]) * (-n[i + 2] + n[i + 4]) * (-n[i + 3] + n[i + 4]) + n[i + 4] * (n[i + 1] * (n[i + 2] - n[i + 4]) * (-n[i + 3] + n[i + 4]) + n[i + 4] * (n[i + 2] * n[i + 3] - (n[i + 2] + n[i + 3]) * n[i + 4]))) * pow(n[i + 5], 3)) / ((n[i + 3] - n[i + 4]) * (-n[i] + n[i + 4]) * (-n[i + 1] + n[i + 4]) * (-n[i + 2] + n[i + 4]) * (n[i + 1] - n[i + 5]) * (n[i + 2] - n[i + 5]) * (n[i + 3] - n[i + 5]));
		weights[s][3][1] = (-4 * n[i] * n[i + 1] * n[i + 2] * n[i + 3] * pow(n[i + 4], 2) + 4 * n[i + 4] * (n[i + 1] * n[i + 2] * n[i + 3] * n[i + 4] + n[i] * (-(n[i + 1] * n[i + 2] * n[i + 3]) + n[i + 2] * n[i + 3] * n[i + 4] + n[i + 1] * (n[i + 2] + n[i + 3]) * n[i + 4])) * n[i + 5] - 4 * (n[i] * n[i + 1] * n[i + 2] * n[i + 3] - (n[i] * n[i + 1] * n[i + 2] + n[i + 1] * n[i + 2] * n[i + 3] + n[i] * (n[i + 1] + n[i + 2]) * n[i + 3]) * n[i + 4] + (n[i + 2] * n[i + 3] + n[i + 1] * (n[i + 2] + n[i + 3]) + n[i] * (n[i + 1] + n[i + 2] + n[i + 3])) * pow(n[i + 4], 2)) * pow(n[i + 5], 2) + 4 * pow(n[i + 4], 3) * pow(n[i + 5], 3)) / ((n[i + 3] - n[i + 4]) * (-n[i] + n[i + 4]) * (-n[i + 1] + n[i + 4]) * (-n[i + 2] + n[i + 4]) * (n[i + 1] - n[i + 5]) * (n[i + 2] - n[i + 5]) * (n[i + 3] - n[i + 5]));
		weights[s][3][2] = (6 * (n[i] * n[i + 1] * n[i + 2] * n[i + 3] * n[i + 4] + (n[i] * n[i + 1] * n[i + 2] * n[i + 3] - (n[i] * n[i + 1] * n[i + 2] + n[i + 1] * n[i + 2] * n[i + 3] + n[i] * (n[i + 1] + n[i + 2]) * n[i + 3]) * n[i + 4]) * n[i + 5] + (n[i] + n[i + 1] + n[i + 2] + n[i + 3] - n[i + 4]) * pow(n[i + 4], 2) * pow(n[i + 5], 2) - pow(n[i + 4], 2) * pow(n[i + 5], 3))) / ((n[i + 3] - n[i + 4]) * (-n[i] + n[i + 4]) * (-n[i + 1] + n[i + 4]) * (-n[i + 2] + n[i + 4]) * (n[i + 1] - n[i + 5]) * (n[i + 2] - n[i + 5]) * (n[i + 3] - n[i + 5]));
		weights[s][3][3] = (-4 * n[i] * n[i + 1] * n[i + 2] * n[i + 3] + 4 * n[i + 4] * (n[i + 1] * n[i + 2] + n[i + 1] * n[i + 3] + n[i + 2] * n[i + 3] + n[i] * (n[i + 1] + n[i + 2] + n[i + 3] - n[i + 4]) - (n[i + 1] + n[i + 2] + n[i + 3]) * n[i + 4] + pow(n[i + 4], 2)) * n[i + 5] - 4 * (n[i] + n[i + 1] + n[i + 2] + n[i + 3] - n[i + 4]) * n[i + 4] * pow(n[i + 5], 2) + 4 * n[i + 4] * pow(n[i + 5], 3)) / ((n[i + 3] - n[i + 4]) * (-n[i] + n[i + 4]) * (-n[i + 1] + n[i + 4]) * (-n[i + 2] + n[i + 4]) * (n[i + 1] - n[i + 5]) * (n[i + 2] - n[i + 5]) * (n[i + 3] - n[i + 5]));
		weights[s][3][4] = 1.0 / ((n[i + 3] - n[i + 4]) * (-n[i] + n[i + 4]) * (-n[i + 1] + n[i + 4]) * (-n[i + 2] + n[i + 4])) + 1 / ((-n[i + 1] + n[i + 2]) * (n[i + 1] - n[i + 3]) * (n[i + 1] - n[i + 4]) * (n[i + 1] - n[i + 5])) + 1 / ((n[i + 1] - n[i + 2]) * (n[i + 2] - n[i + 3]) * (n[i + 2] - n[i + 4]) * (n[i + 2] - n[i + 5])) + 1 / ((n[i + 1] - n[i + 3]) * (-n[i + 2] + n[i + 3]) * (n[i + 3] - n[i + 4]) * (n[i + 3] - n[i + 5]));
		//Part 5:
		weights[s][4][0] = pow(n[i + 5], 4) / ((-n[i + 1] + n[i + 5]) * (-n[i + 2] + n[i + 5]) * (-n[i + 3] + n[i + 5]) * (-n[i + 4] + n[i + 5]));
		weights[s][4][1] = (4 * pow(n[i + 5], 3)) / ((n[i + 4] - n[i + 5]) * (-n[i + 1] + n[i + 5]) * (-n[i + 2] + n[i + 5]) * (-n[i + 3] + n[i + 5]));
		weights[s][4][2] = (6 * pow(n[i + 5], 2)) / ((-n[i + 1] + n[i + 5]) * (-n[i + 2] + n[i + 5]) * (-n[i + 3] + n[i + 5]) * (-n[i + 4] + n[i + 5]));
		weights[s][4][3] = (4 * n[i + 5]) / ((n[i + 4] - n[i + 5]) * (-n[i + 1] + n[i + 5]) * (-n[i + 2] + n[i + 5]) * (-n[i + 3] + n[i + 5]));
		weights[s][4][4] = 1.0 / ((n[i + 1] - n[i + 5]) * (n[i + 2] - n[i + 5]) * (n[i + 3] - n[i + 5]) * (n[i + 4] - n[i + 5]));
	}
	return weights;
}

void SetBoundaryConditions0_1D_OR_1(vector<double> n, vector<vector<double> >& bcs)
{
	bcs.resize(0);
}

void SetBoundaryConditions0_1D_CO_1(vector<double> n, vector<vector<double> >& bcs, double rc)
{
	bcs.resize(0);
}

void SetBoundaryConditions1_1(vector<double> n, vector<vector<double> >& bcs)
{
	vector<double> bc(4);

	int i = n.size() - 6;

	bc[0] = 1.0;
	bc[1] = (-n[i + 3] + n[i + 4]) / (n[i + 1] - n[i + 3]);
	bc[2] = ((n[i + 3] - n[i + 4]) * (-n[i + 3] + n[i + 5])) / ((n[i + 2] - n[i + 3]) * (-n[i + 1] + n[i + 3]));
	bc[3] = 0.0;
	bcs.push_back(bc);
}

void SetBoundaryConditions1_2(vector<double> n, vector<vector<double> >& bcs)
{
	vector<double> bc(4);

	int i = n.size() - 6;

	bc[0] = 1.0;
	bc[1] = 0.0;
	bc[2] = ((n[i + 3] - n[i + 4]) * (n[i + 2] - n[i + 5])) / ((n[i + 2] - n[i + 3]) * (n[i + 1] - n[i + 4]));
	bc[3] = 0.0;
	bcs.push_back(bc);

	bc[0] = 0.0;
	bc[1] = 1.0;
	bc[2] = (n[i + 1] * (n[i + 2] - n[i + 3]) - n[i + 2] * n[i + 3] + n[i + 3] * n[i + 4] + n[i + 3] * n[i + 5] - n[i + 4] * n[i + 5]) / ((n[i + 2] - n[i + 3]) * (n[i + 1] - n[i + 4]));
	bc[3] = 0.0;
	bcs.push_back(bc);
}

void SetBoundaryConditions2_1(vector<double> n, vector<vector<double> >& bcs)
{
	vector<double> bc(4);

	int i = n.size() - 6;

	bc[0] = 1.0;
	bc[1] = 1.0;
	bc[2] = 1.0;
	bc[3] = 0.0;
	bcs.push_back(bc);
}

void SetBoundaryConditions2_2(vector<double> n, vector<vector<double> >& bcs)
{
	vector<double> bc(4);

	int i = n.size() - 6;

	bc[0] = 1.0;
	bc[1] = 0.0;
	bc[2] = (-n[i + 2] + n[i + 3]) / (n[i + 1] - n[i + 4]);
	bc[3] = 0.0;
	bcs.push_back(bc);

	bc[0] = 0.0;
	bc[1] = 1.0;
	bc[2] = (n[i + 2] - n[i + 4]) * (1 / (n[i + 1] - n[i + 3]) + 1 / (n[i + 2] - n[i + 5]));
	bc[3] = 0.0;
	bcs.push_back(bc);
}

void SetBoundaryConditions1_MM_1(vector<double> n, vector<vector<double> >& bcs, double rc, double paramMcMillanExponent)
{
	vector<double> bc(2);

	int i = -1;

	//mcMillanParameter: spline0, spline1
	bc[0] = (pow(rc, -1 + paramMcMillanExponent) * ((pow(n[i + 3], 2) * (-n[i + 2] + n[i + 4]) + (pow(n[i + 3], 2) + n[i + 2] * n[i + 4] - 2 * n[i + 3] * n[i + 4]) * n[i + 5] - n[i + 1] * (pow(n[i + 3], 2) - n[i + 4] * n[i + 5] + n[i + 2] * (-2 * n[i + 3] + n[i + 4] + n[i + 5]))) * paramMcMillanExponent + 3 * (n[i + 2] * n[i + 3] + n[i + 1] * (-n[i + 2] + n[i + 3]) - n[i + 3] * n[i + 4] - n[i + 3] * n[i + 5] + n[i + 4] * n[i + 5]) * rc)) / (3. * (n[i + 3] - n[i + 4]) * (n[i + 3] - n[i + 5]));
	bc[1] = ((-n[i + 2] + n[i + 5]) * ((n[i + 3] - n[i + 4]) * paramMcMillanExponent - 3 * rc) * pow(rc, -1 + paramMcMillanExponent)) / (3. * (n[i + 3] - n[i + 5]));
	bcs.push_back(bc);

	//spline2Parameter: spline0, spline1
	bc[0] = ((n[i + 1] - n[i + 3]) * (n[i + 2] - n[i + 3])) / ((n[i + 3] - n[i + 4]) * (n[i + 3] - n[i + 5]));
	bc[1] = (-n[i + 2] + n[i + 3]) / (n[i + 3] - n[i + 5]);
	bcs.push_back(bc);
}

void SetBoundaryConditions1_EXP_1(vector<double> n, vector<vector<double> >& bcs, double rc)
{
	//TODO: recalculate expressions in mathematica notebook!
	vector<double> bc(4);

	int i = n.size() - 6;

	bc[0] = ((n[i + 2] - n[i + 3]) * (-n[i + 1] + n[i + 3]) * (n[i + 3] - n[i + 4]) * (n[i + 2] - n[i + 5]) * (n[i + 2] + n[i + 3] - n[i + 4] - n[i + 5])) / ((n[i + 1] - n[i + 4]) * (n[i + 2] - n[i + 4]) * (-((n[i + 2] - n[i + 3]) * n[i + 4] * (n[i + 2] + n[i + 4])) + n[i + 4] * (2 * n[i + 2] - 3 * n[i + 3] + n[i + 4]) * n[i + 5] + (n[i + 3] - n[i + 4]) * pow(n[i + 5], 2) + n[i + 1] * (pow(n[i + 2], 2) - pow(n[i + 4], 2) + n[i + 2] * (-n[i + 3] + n[i + 4] - 2 * n[i + 5]) + (n[i + 3] + n[i + 4]) * n[i + 5])));
	bc[1] = ((n[i + 1] - n[i + 3]) * (n[i + 2] - n[i + 3]) * (n[i + 2] - n[i + 5])) / (-((n[i + 2] - n[i + 3]) * n[i + 4] * (n[i + 2] + n[i + 4])) + n[i + 4] * (2 * n[i + 2] - 3 * n[i + 3] + n[i + 4]) * n[i + 5] + (n[i + 3] - n[i + 4]) * pow(n[i + 5], 2) + n[i + 1] * (pow(n[i + 2], 2) - pow(n[i + 4], 2) + n[i + 2] * (-n[i + 3] + n[i + 4] - 2 * n[i + 5]) + (n[i + 3] + n[i + 4]) * n[i + 5]));
	bc[2] = ((n[i + 1] - n[i + 3]) * (n[i + 2] - n[i + 3]) * (n[i + 2] - n[i + 5]) * (n[i + 2] - n[i + 5] + 3 * rc)) / (3. * (-((n[i + 2] - n[i + 3]) * n[i + 4] * (n[i + 2] + n[i + 4])) + n[i + 4] * (2 * n[i + 2] - 3 * n[i + 3] + n[i + 4]) * n[i + 5] + (n[i + 3] - n[i + 4]) * pow(n[i + 5], 2) + n[i + 1] * (pow(n[i + 2], 2) - pow(n[i + 4], 2) + n[i + 2] * (-n[i + 3] + n[i + 4] - 2 * n[i + 5]) + (n[i + 3] + n[i + 4]) * n[i + 5])));
	bc[3] = ((n[i + 1] - n[i + 3]) * (n[i + 2] - n[i + 3]) * (n[i + 2] - n[i + 5]) * (n[i + 2] - n[i + 5] + 3 * rc * log(rc))) / (3. * (-((n[i + 2] - n[i + 3]) * n[i + 4] * (n[i + 2] + n[i + 4])) + n[i + 4] * (2 * n[i + 2] - 3 * n[i + 3] + n[i + 4]) * n[i + 5] + (n[i + 3] - n[i + 4]) * pow(n[i + 5], 2) + n[i + 1] * (pow(n[i + 2], 2) - pow(n[i + 4], 2) + n[i + 2] * (-n[i + 3] + n[i + 4] - 2 * n[i + 5]) + (n[i + 3] + n[i + 4]) * n[i + 5])) * rc);
	bcs.push_back(bc);

	bc[0] = ((n[i + 2] - n[i + 3]) * (n[i + 3] - n[i + 4]) * ((n[i + 2] - n[i + 3]) * (n[i + 3] - 2 * n[i + 4]) * n[i + 4] + (n[i + 3] * (n[i + 2] + n[i + 3]) - 4 * n[i + 3] * n[i + 4] + 2 * pow(n[i + 4], 2)) * n[i + 5] + n[i + 1] * ((n[i + 3] - 2 * n[i + 4]) * n[i + 4] + n[i + 3] * n[i + 5] - n[i + 2] * (2 * n[i + 3] - 3 * n[i + 4] + n[i + 5])))) / ((n[i + 1] - n[i + 4]) * (n[i + 2] - n[i + 5]) * (-((n[i + 2] - n[i + 3]) * n[i + 4] * (n[i + 2] + n[i + 4])) + n[i + 4] * (2 * n[i + 2] - 3 * n[i + 3] + n[i + 4]) * n[i + 5] + (n[i + 3] - n[i + 4]) * pow(n[i + 5], 2) + n[i + 1] * (pow(n[i + 2], 2) - pow(n[i + 4], 2) + n[i + 2] * (-n[i + 3] + n[i + 4] - 2 * n[i + 5]) + (n[i + 3] + n[i + 4]) * n[i + 5])));
	bc[1] = ((n[i + 2] - n[i + 3]) * (n[i + 2] - n[i + 4]) * (n[i + 1] * (n[i + 2] - n[i + 4]) - n[i + 2] * n[i + 4] + n[i + 3] * n[i + 4] - n[i + 3] * n[i + 5] + n[i + 4] * n[i + 5])) / ((n[i + 2] - n[i + 5]) * (-((n[i + 2] - n[i + 3]) * n[i + 4] * (n[i + 2] + n[i + 4])) + n[i + 4] * (2 * n[i + 2] - 3 * n[i + 3] + n[i + 4]) * n[i + 5] + (n[i + 3] - n[i + 4]) * pow(n[i + 5], 2) + n[i + 1] * (pow(n[i + 2], 2) - pow(n[i + 4], 2) + n[i + 2] * (-n[i + 3] + n[i + 4] - 2 * n[i + 5]) + (n[i + 3] + n[i + 4]) * n[i + 5])));
	bc[2] = ((n[i + 2] - n[i + 3]) * (n[i + 2] - n[i + 4]) * (-(n[i + 4] * (n[i + 3] * n[i + 4] - 2 * n[i + 3] * n[i + 5] + n[i + 4] * n[i + 5])) + 3 * (n[i + 3] * n[i + 4] - n[i + 3] * n[i + 5] + n[i + 4] * n[i + 5]) * rc + n[i + 2] * (pow(n[i + 4], 2) - n[i + 3] * n[i + 5] - 3 * n[i + 4] * rc) + n[i + 1] * (pow(n[i + 4], 2) - n[i + 3] * n[i + 5] - 3 * n[i + 4] * rc + n[i + 2] * (n[i + 3] - 2 * n[i + 4] + n[i + 5] + 3 * rc)))) / (3. * (n[i + 2] - n[i + 5]) * (-((n[i + 2] - n[i + 3]) * n[i + 4] * (n[i + 2] + n[i + 4])) + n[i + 4] * (2 * n[i + 2] - 3 * n[i + 3] + n[i + 4]) * n[i + 5] + (n[i + 3] - n[i + 4]) * pow(n[i + 5], 2) + n[i + 1] * (pow(n[i + 2], 2) - pow(n[i + 4], 2) + n[i + 2] * (-n[i + 3] + n[i + 4] - 2 * n[i + 5]) + (n[i + 3] + n[i + 4]) * n[i + 5])));
	bc[3] = ((n[i + 2] - n[i + 3]) * (n[i + 2] - n[i + 4]) * ((n[i + 2] - n[i + 3]) * pow(n[i + 4], 2) - (n[i + 2] * n[i + 3] + n[i + 4] * (-2 * n[i + 3] + n[i + 4])) * n[i + 5] + n[i + 1] * (pow(n[i + 4], 2) - n[i + 3] * n[i + 5] + n[i + 2] * (n[i + 3] - 2 * n[i + 4] + n[i + 5])) + 3 * (n[i + 1] * (n[i + 2] - n[i + 4]) - n[i + 2] * n[i + 4] + n[i + 3] * n[i + 4] - n[i + 3] * n[i + 5] + n[i + 4] * n[i + 5]) * rc * log(rc))) / (3. * (n[i + 2] - n[i + 5]) * (-((n[i + 2] - n[i + 3]) * n[i + 4] * (n[i + 2] + n[i + 4])) + n[i + 4] * (2 * n[i + 2] - 3 * n[i + 3] + n[i + 4]) * n[i + 5] + (n[i + 3] - n[i + 4]) * pow(n[i + 5], 2) + n[i + 1] * (pow(n[i + 2], 2) - pow(n[i + 4], 2) + n[i + 2] * (-n[i + 3] + n[i + 4] - 2 * n[i + 5]) + (n[i + 3] + n[i + 4]) * n[i + 5])) * rc);
	bcs.push_back(bc);
}

void SetBoundaryConditions1_EXP_2(vector<double> n, vector<vector<double> >& bcs, double rc)
{
	vector<double> bc(2);

	int i = n.size() - 6;

	//secondSecondLastSpline: secondLastSpline, lastSpline
	bc[0] = (-n[i + 3] + n[i + 4]) / (n[i + 1] - n[i + 3]);
	bc[1] = -(((n[i + 3] - n[i + 4]) * (n[i + 3] - n[i + 5])) / ((n[i + 2] - n[i + 3]) * (-n[i + 1] + n[i + 3])));
	bcs.push_back(bc);

	//constantParameter: secondLastSpline, lastSpline
	bc[0] = (n[i + 1] - n[i + 4]) / (n[i + 1] - n[i + 3]);
	bc[1] = (n[i + 2] * n[i + 3] + n[i + 1] * (-n[i + 2] + n[i + 3]) - n[i + 3] * n[i + 4] - n[i + 3] * n[i + 5] + n[i + 4] * n[i + 5]) / ((n[i + 2] - n[i + 3]) * (-n[i + 1] + n[i + 3]));
	bcs.push_back(bc);

	//linearParameter: secondLastSpline, lastSpline
	bc[0] = ((n[i + 1] - n[i + 4]) * (n[i + 2] - n[i + 3] + 3 * rc)) / (3. * (n[i + 1] - n[i + 3]));
	bc[1] = (n[i + 3] * (-2 * n[i + 4] * n[i + 5] + n[i + 3] * (n[i + 4] + n[i + 5])) + 3 * n[i + 4] * n[i + 5] * rc - 3 * n[i + 3] * (n[i + 4] + n[i + 5]) * rc + n[i + 2] * (-pow(n[i + 3], 2) + n[i + 4] * n[i + 5] + 3 * n[i + 3] * rc) - n[i + 1] * (pow(n[i + 3], 2) - n[i + 4] * n[i + 5] - 3 * n[i + 3] * rc + n[i + 2] * (-2 * n[i + 3] + n[i + 4] + n[i + 5] + 3 * rc))) / (3. * (n[i + 2] - n[i + 3]) * (-n[i + 1] + n[i + 3]));
	bcs.push_back(bc);
}

void SetBoundaryConditions1_EXP_3(vector<double> n, vector<vector<double> >& bcs, double rc)
{
	vector<double> bc(2);

	int i = n.size() - 6;

	//secondSecondLastSpline: secondLastSpline, lastSpline
	bc[0] = (-n[i + 3] + n[i + 4]) / (n[i + 1] - n[i + 3]);
	bc[1] = -(((n[i + 3] - n[i + 4]) * (n[i + 3] - n[i + 5])) / ((n[i + 2] - n[i + 3]) * (-n[i + 1] + n[i + 3])));
	bcs.push_back(bc);

	//constantParameter: secondLastSpline, lastSpline
	bc[0] = (n[i + 1] - n[i + 4]) / (n[i + 1] - n[i + 3]);
	bc[1] = (n[i + 2] * n[i + 3] + n[i + 1] * (-n[i + 2] + n[i + 3]) - n[i + 3] * n[i + 4] - n[i + 3] * n[i + 5] + n[i + 4] * n[i + 5]) / ((n[i + 2] - n[i + 3]) * (-n[i + 1] + n[i + 3]));
	bcs.push_back(bc);

	//linearParameter: secondLastSpline, lastSpline
	bc[0] = ((n[i + 1] - n[i + 4]) * (n[i + 2] - n[i + 3] + 3 * rc)) / (3. * (n[i + 1] - n[i + 3]));
	bc[1] = (n[i + 3] * (-2 * n[i + 4] * n[i + 5] + n[i + 3] * (n[i + 4] + n[i + 5])) + 3 * n[i + 4] * n[i + 5] * rc - 3 * n[i + 3] * (n[i + 4] + n[i + 5]) * rc + n[i + 2] * (-pow(n[i + 3], 2) + n[i + 4] * n[i + 5] + 3 * n[i + 3] * rc) - n[i + 1] * (pow(n[i + 3], 2) - n[i + 4] * n[i + 5] - 3 * n[i + 3] * rc + n[i + 2] * (-2 * n[i + 3] + n[i + 4] + n[i + 5] + 3 * rc))) / (3. * (n[i + 2] - n[i + 3]) * (-n[i + 1] + n[i + 3]));
	bcs.push_back(bc);

	//logParameter: secondLastSpline, lastSpline
	bc[0] = ((n[i + 1] - n[i + 4]) * (n[i + 2] - n[i + 3] + 3 * rc * log(rc))) / (3. * (n[i + 1] - n[i + 3]) * rc);
	bc[1] = (pow(n[i + 3], 2) * (-n[i + 2] + n[i + 4]) + (pow(n[i + 3], 2) + n[i + 2] * n[i + 4] - 2 * n[i + 3] * n[i + 4]) * n[i + 5] - n[i + 1] * (pow(n[i + 3], 2) - n[i + 4] * n[i + 5] + n[i + 2] * (-2 * n[i + 3] + n[i + 4] + n[i + 5])) + 3 * (n[i + 2] * n[i + 3] + n[i + 1] * (-n[i + 2] + n[i + 3]) - n[i + 3] * n[i + 4] - n[i + 3] * n[i + 5] + n[i + 4] * n[i + 5]) * rc * log(rc)) / (3. * (n[i + 2] - n[i + 3]) * (-n[i + 1] + n[i + 3]) * rc);
	bcs.push_back(bc);
}

void SetBoundaryConditions1_MM_1_4thorder(vector<double> n, vector<vector<double> >& bcs, double rc, double paramMcMillanExponent)
{
	vector<double> bc(3);

	int i = -1;

	//mcMillanParameter: spline0, spline1, spline 2
	bc[0] = (pow(rc, -2 + paramMcMillanExponent)
			* ((pow(n[i + 4], 2) * (2 * n[i + 4] * n[i + 5] * n[i + 6] - pow(n[i + 4], 2) * (n[i + 5] + n[i + 6]) + n[i + 3] * (pow(n[i + 4], 2) - n[i + 5] * n[i + 6]) + n[i + 2] * (pow(n[i + 4], 2) - n[i + 5] * n[i + 6] + n[i + 3] * (-2 * n[i + 4] + n[i + 5] + n[i + 6]))) - (pow(n[i + 4], 2) * (-(n[i + 2] * n[i + 3]) + pow(n[i + 4], 2) + (n[i + 2] + n[i + 3] - 2 * n[i + 4]) * n[i + 5]) + ((n[i + 2] + n[i + 3] - 2 * n[i + 4]) * pow(n[i + 4], 2) + (n[i + 2] * n[i + 3] - 2 * (n[i + 2] + n[i + 3]) * n[i + 4] + 3 * pow(n[i + 4], 2)) * n[i + 5]) * n[i + 6]) * n[i + 7]
					+ n[i + 1] * (pow(n[i + 4], 2) * (pow(n[i + 4], 2) - n[i + 5] * n[i + 6] + n[i + 3] * (-2 * n[i + 4] + n[i + 5] + n[i + 6])) + pow(n[i + 4], 2) * (n[i + 3] - n[i + 5]) * n[i + 7] - (pow(n[i + 4], 2) + n[i + 3] * n[i + 5] - 2 * n[i + 4] * n[i + 5]) * n[i + 6] * n[i + 7] + n[i + 2] * (-2 * pow(n[i + 4], 3) - n[i + 5] * n[i + 6] * n[i + 7] + pow(n[i + 4], 2) * (n[i + 5] + n[i + 6] + n[i + 7]) + n[i + 3] * (3 * pow(n[i + 4], 2) + n[i + 5] * n[i + 6] + (n[i + 5] + n[i + 6]) * n[i + 7] - 2 * n[i + 4] * (n[i + 5] + n[i + 6] + n[i + 7]))))) * pow(paramMcMillanExponent, 2)
					+ 12 * (n[i + 1] * (n[i + 2] - n[i + 4]) * (n[i + 3] - n[i + 4]) + n[i + 4] * (n[i + 3] * n[i + 4] + n[i + 2] * (-n[i + 3] + n[i + 4]) - n[i + 4] * n[i + 5] - n[i + 4] * n[i + 6] + n[i + 5] * n[i + 6]) - (n[i + 4] - n[i + 5]) * (n[i + 4] - n[i + 6]) * n[i + 7]) * pow(rc, 2)
					+ paramMcMillanExponent
							* (2 * n[i + 2] * n[i + 3] * pow(n[i + 4], 3) - n[i + 2] * pow(n[i + 4], 4) - n[i + 3] * pow(n[i + 4], 4) - n[i + 2] * n[i + 3] * pow(n[i + 4], 2) * n[i + 5] + pow(n[i + 4], 4) * n[i + 5] - n[i + 2] * n[i + 3] * pow(n[i + 4], 2) * n[i + 6] + pow(n[i + 4], 4) * n[i + 6] + n[i + 2] * pow(n[i + 4], 2) * n[i + 5] * n[i + 6] + n[i + 3] * pow(n[i + 4], 2) * n[i + 5] * n[i + 6] - 2 * pow(n[i + 4], 3) * n[i + 5] * n[i + 6] - n[i + 2] * n[i + 3] * pow(n[i + 4], 2) * n[i + 7] + pow(n[i + 4], 4) * n[i + 7] + n[i + 2] * pow(n[i + 4], 2) * n[i + 5] * n[i + 7] + n[i + 3] * pow(n[i + 4], 2) * n[i + 5] * n[i + 7] - 2 * pow(n[i + 4], 3) * n[i + 5] * n[i + 7] + n[i + 2] * pow(n[i + 4], 2) * n[i + 6] * n[i + 7] + n[i + 3] * pow(n[i + 4], 2) * n[i + 6] * n[i + 7]
									- 2 * pow(n[i + 4], 3) * n[i + 6] * n[i + 7] + n[i + 2] * n[i + 3] * n[i + 5] * n[i + 6] * n[i + 7] - 2 * n[i + 2] * n[i + 4] * n[i + 5] * n[i + 6] * n[i + 7] - 2 * n[i + 3] * n[i + 4] * n[i + 5] * n[i + 6] * n[i + 7] + 3 * pow(n[i + 4], 2) * n[i + 5] * n[i + 6] * n[i + 7]
									+ 3 * (n[i + 4] * (-2 * n[i + 2] * pow(n[i + 4], 2) + n[i + 2] * n[i + 5] * n[i + 6] - 3 * n[i + 4] * n[i + 5] * n[i + 6] + 2 * pow(n[i + 4], 2) * (n[i + 5] + n[i + 6]) - n[i + 2] * n[i + 3] * (-3 * n[i + 4] + n[i + 5] + n[i + 6]) + n[i + 3] * (-2 * pow(n[i + 4], 2) + n[i + 5] * n[i + 6])) + (-(n[i + 2] * n[i + 3] * n[i + 4]) + 2 * pow(n[i + 4], 3) + (n[i + 2] + n[i + 3] - 3 * n[i + 4]) * n[i + 4] * n[i + 5] + (n[i + 2] + n[i + 3] - 3 * n[i + 4]) * (n[i + 4] - n[i + 5]) * n[i + 6]) * n[i + 7]) * rc
									- n[i + 1]
											* (-2 * n[i + 3] * pow(n[i + 4], 3) + pow(n[i + 4], 4) + n[i + 3] * pow(n[i + 4], 2) * n[i + 5] + n[i + 3] * pow(n[i + 4], 2) * n[i + 6] - pow(n[i + 4], 2) * n[i + 5] * n[i + 6] + n[i + 3] * pow(n[i + 4], 2) * n[i + 7] - pow(n[i + 4], 2) * n[i + 5] * n[i + 7] - pow(n[i + 4], 2) * n[i + 6] * n[i + 7] - n[i + 3] * n[i + 5] * n[i + 6] * n[i + 7] + 2 * n[i + 4] * n[i + 5] * n[i + 6] * n[i + 7] + 3 * (2 * pow(n[i + 4], 3) + n[i + 5] * n[i + 6] * n[i + 7] + n[i + 3] * n[i + 4] * (-3 * n[i + 4] + n[i + 5] + n[i + 6] + n[i + 7]) - n[i + 4] * (n[i + 5] * n[i + 6] + (n[i + 5] + n[i + 6]) * n[i + 7])) * rc
													+ n[i + 2] * (-2 * pow(n[i + 4], 3) - n[i + 5] * n[i + 6] * n[i + 7] + pow(n[i + 4], 2) * (n[i + 5] + n[i + 6] + n[i + 7] - 9 * rc) + 3 * n[i + 4] * (n[i + 5] + n[i + 6] + n[i + 7]) * rc + n[i + 3] * (3 * pow(n[i + 4], 2) + n[i + 5] * n[i + 6] + n[i + 5] * n[i + 7] + n[i + 6] * n[i + 7] - 2 * n[i + 4] * (n[i + 5] + n[i + 6] + n[i + 7]) + 9 * n[i + 4] * rc - 3 * (n[i + 5] + n[i + 6] + n[i + 7]) * rc)))))) / (12. * (n[i + 4] - n[i + 5]) * (n[i + 4] - n[i + 6]) * (n[i + 4] - n[i + 7]));
	bc[1] = (pow(rc, -2 + paramMcMillanExponent)
			* ((n[i + 4] - n[i + 5]) * (2 * n[i + 4] * n[i + 6] * n[i + 7] - pow(n[i + 4], 2) * (n[i + 6] + n[i + 7]) + n[i + 3] * (pow(n[i + 4], 2) - n[i + 6] * n[i + 7]) + n[i + 2] * (pow(n[i + 4], 2) - n[i + 6] * n[i + 7] + n[i + 3] * (-2 * n[i + 4] + n[i + 6] + n[i + 7]))) * (-1 + paramMcMillanExponent) * paramMcMillanExponent + 3 * (-(n[i + 4] * (2 * n[i + 4] - n[i + 5]) * (n[i + 3] - n[i + 6])) + (2 * pow(n[i + 4], 2) + (n[i + 3] + n[i + 5]) * n[i + 6] - n[i + 4] * (n[i + 5] + 3 * n[i + 6])) * n[i + 7] + n[i + 2] * (n[i + 4] * (-2 * n[i + 4] + n[i + 5]) + n[i + 3] * (3 * n[i + 4] - n[i + 5] - n[i + 6] - n[i + 7]) + n[i + 6] * n[i + 7])) * paramMcMillanExponent * rc
					+ 12 * (n[i + 3] * n[i + 4] + n[i + 2] * (-n[i + 3] + n[i + 4]) - n[i + 4] * n[i + 6] - n[i + 4] * n[i + 7] + n[i + 6] * n[i + 7]) * pow(rc, 2))) / (12. * (n[i + 4] - n[i + 6]) * (n[i + 4] - n[i + 7]));
	bc[2] = ((n[i + 3] - n[i + 7]) * pow(rc, -2 + paramMcMillanExponent) * (-((n[i + 4] - n[i + 5]) * (n[i + 4] - n[i + 6]) * (-1 + paramMcMillanExponent) * paramMcMillanExponent) + 3 * (2 * n[i + 4] - n[i + 5] - n[i + 6]) * paramMcMillanExponent * rc - 12 * pow(rc, 2))) / (12. * (-n[i + 4] + n[i + 7]));
	bcs.push_back(bc);

	//spline2Parameter: spline0, spline1, spline 2
	bc[0] = ((n[i + 1] - n[i + 4]) * (-n[i + 2] + n[i + 4]) * (-n[i + 3] + n[i + 4])) / ((n[i + 4] - n[i + 5]) * (n[i + 4] - n[i + 6]) * (-n[i + 4] + n[i + 7]));
	bc[1] = ((n[i + 2] - n[i + 4]) * (n[i + 3] - n[i + 4])) / ((n[i + 4] - n[i + 6]) * (n[i + 4] - n[i + 7]));
	bc[2] = (-n[i + 3] + n[i + 4]) / (n[i + 4] - n[i + 7]);
	bcs.push_back(bc);
}

void SetBoundaryConditions1_EXP_2_4thorder(vector<double> n, vector<vector<double> >& bcs, double rc)
{
	vector<double> bc(3);

	int i = n.size() - 8;

	//-4th spline: -3rd spline, -2nd splin, last spline
	bc[0] = (-n[i + 4] + n[i + 5]) / (n[i + 1] - n[i + 4]);
	bc[1] = -(((n[i + 4] - n[i + 5]) * (n[i + 4] - n[i + 6])) / ((n[i + 2] - n[i + 4]) * (-n[i + 1] + n[i + 4])));
	bc[2] = -(((n[i + 4] - n[i + 5]) * (n[i + 4] - n[i + 6]) * (n[i + 4] - n[i + 7])) / ((n[i + 3] - n[i + 4]) * (-n[i + 1] + n[i + 4]) * (-n[i + 2] + n[i + 4])));
	bcs.push_back(bc);

	//constantParameter: -3rd spline, -2nd splin, last spline
	bc[0] = (n[i + 1] - n[i + 5]) / (n[i + 1] - n[i + 4]);
	bc[1] = (n[i + 2] * n[i + 4] + n[i + 1] * (-n[i + 2] + n[i + 4]) - n[i + 4] * n[i + 5] - n[i + 4] * n[i + 6] + n[i + 5] * n[i + 6]) / ((n[i + 1] - n[i + 4]) * (-n[i + 2] + n[i + 4]));
	bc[2] = (n[i + 1] * (n[i + 2] - n[i + 4]) * (n[i + 3] - n[i + 4]) + n[i + 4] * (n[i + 3] * n[i + 4] + n[i + 2] * (-n[i + 3] + n[i + 4]) - n[i + 4] * n[i + 5] - n[i + 4] * n[i + 6] + n[i + 5] * n[i + 6]) - (n[i + 4] - n[i + 5]) * (n[i + 4] - n[i + 6]) * n[i + 7]) / ((n[i + 1] - n[i + 4]) * (-n[i + 2] + n[i + 4]) * (-n[i + 3] + n[i + 4]));
	bcs.push_back(bc);

	//linearParameter: -3rd spline, -2nd splin, last spline
	bc[0] = ((n[i + 1] - n[i + 5]) * (n[i + 2] + n[i + 3] - 2 * n[i + 4] + 4 * rc)) / (4. * (n[i + 1] - n[i + 4]));
	bc[1] = -(n[i + 3] * n[i + 4] * n[i + 5] - 2 * pow(n[i + 4], 2) * n[i + 5] + n[i + 3] * n[i + 4] * n[i + 6] - 2 * pow(n[i + 4], 2) * n[i + 6] - n[i + 3] * n[i + 5] * n[i + 6] + 3 * n[i + 4] * n[i + 5] * n[i + 6] + 4 * (-(n[i + 5] * n[i + 6]) + n[i + 4] * (n[i + 5] + n[i + 6])) * rc - n[i + 2] * (n[i + 5] * n[i + 6] + n[i + 4] * (n[i + 3] - 2 * n[i + 4] + 4 * rc)) + n[i + 1] * (-(n[i + 5] * n[i + 6]) + n[i + 4] * (-n[i + 3] + 2 * n[i + 4] - 4 * rc) + n[i + 2] * (n[i + 3] - 3 * n[i + 4] + n[i + 5] + n[i + 6] + 4 * rc))) / (4. * (n[i + 1] - n[i + 4]) * (-n[i + 2] + n[i + 4]));
	bc[2] = (-2 * n[i + 3] * pow(n[i + 4], 3) + 2 * pow(n[i + 4], 3) * n[i + 5] + 2 * pow(n[i + 4], 3) * n[i + 6] + n[i + 3] * n[i + 4] * n[i + 5] * n[i + 6] - 3 * pow(n[i + 4], 2) * n[i + 5] * n[i + 6] + 2 * pow(n[i + 4], 3) * n[i + 7] + n[i + 3] * n[i + 4] * n[i + 5] * n[i + 7] - 3 * pow(n[i + 4], 2) * n[i + 5] * n[i + 7] + n[i + 3] * n[i + 4] * n[i + 6] * n[i + 7] - 3 * pow(n[i + 4], 2) * n[i + 6] * n[i + 7] - n[i + 3] * n[i + 5] * n[i + 6] * n[i + 7] + 3 * n[i + 4] * n[i + 5] * n[i + 6] * n[i + 7] + 4 * (n[i + 3] * pow(n[i + 4], 2) + n[i + 4] * n[i + 5] * n[i + 6] - n[i + 5] * n[i + 6] * n[i + 7] + n[i + 4] * (n[i + 5] + n[i + 6]) * n[i + 7] - pow(n[i + 4], 2) * (n[i + 5] + n[i + 6] + n[i + 7])) * rc
			+ n[i + 2] * (-2 * pow(n[i + 4], 3) + n[i + 4] * n[i + 5] * n[i + 6] - n[i + 5] * n[i + 6] * n[i + 7] + n[i + 4] * (n[i + 5] + n[i + 6]) * n[i + 7] + n[i + 3] * n[i + 4] * (3 * n[i + 4] - n[i + 5] - n[i + 6] - n[i + 7] - 4 * rc) + 4 * pow(n[i + 4], 2) * rc) + n[i + 1] * (-2 * pow(n[i + 4], 3) + n[i + 4] * n[i + 5] * n[i + 6] + n[i + 4] * n[i + 5] * n[i + 7] + n[i + 4] * n[i + 6] * n[i + 7] - n[i + 5] * n[i + 6] * n[i + 7] - n[i + 2] * (n[i + 3] - n[i + 4]) * (3 * n[i + 4] - n[i + 5] - n[i + 6] - n[i + 7] - 4 * rc) + n[i + 3] * n[i + 4] * (3 * n[i + 4] - n[i + 5] - n[i + 6] - n[i + 7] - 4 * rc) + 4 * pow(n[i + 4], 2) * rc)) / (4. * (n[i + 1] - n[i + 4]) * (-n[i + 2] + n[i + 4]) * (-n[i + 3] + n[i + 4]));
	bcs.push_back(bc);
}

void SetBoundaryConditions3_1D_OR_1(vector<double> n, vector<vector<double> >& bcs, bool uniform)
{
	//no bc
	vector<double> bc(3);

	if (uniform)
	{
		bc[0] = 1.0;
		bc[1] = 0.0;
		bc[2] = 0.0;
		bcs.push_back(bc);

		bc[0] = 0.0;
		bc[1] = 1.0;
		bc[2] = 0.0;
		bcs.push_back(bc);

		bc[0] = 0.0;
		bc[1] = 0.0;
		bc[2] = 1.0;
		bcs.push_back(bc);
	}
	else
	{
		int i = n.size() - 6;

		//first param: 1st spline, 2nd spline, 3rd spline
		bc[0] = 1.0;
		bc[1] = 0.0;
		bc[2] = 0.0;
		bcs.push_back(bc);

		//second param: 1st spline, 2nd spline, 3rd spline
		bc[0] = 0.0;
		bc[1] = 1.0;
		bc[2] = 0.0;
		bcs.push_back(bc);

		//thirdparam: 1st spline, 2nd spline, 3rd spline
		bc[0] = 0.0;
		bc[1] = 0.0;
		bc[2] = 1.0;
		bcs.push_back(bc);
	}
}

void SetBoundaryConditions3_1D_OR_2(vector<double> n, vector<vector<double> >& bcs, bool uniform)
{
	//first derivatie = 0
	vector<double> bc(3);

	if (uniform)
	{
		bc[0] = 0.0;
		bc[1] = 1.0;
		bc[2] = 0.0;
		bcs.push_back(bc);

		bc[0] = 1.0;
		bc[1] = 0.0;
		bc[2] = 1.0;
		bcs.push_back(bc);
	}
	else
	{
		int i = n.size() - 6;

		//first param: 1st spline, 2nd spline, 3rd spline
		bc[0] = 0.0;
		bc[1] = 1.0;
		bc[2] = 0.0;
		bcs.push_back(bc);

		//second param: 1st spline, 2nd spline, 3rd spline
		bc[0] = 1.0;
		bc[1] = 0.0;
		bc[2] = 1.0;
		bcs.push_back(bc);
	}
}

void SetBoundaryConditions3_1D_OR_3(vector<double> n, vector<vector<double> >& bcs, bool uniform)
{
	//second derivative = 0
	vector<double> bc(3);

	if (uniform)
	{
		bc[0] = -1.0;
		bc[1] = 0.0;
		bc[2] = 1.0;
		bcs.push_back(bc);

		bc[0] = 2.0;
		bc[1] = 1.0;
		bc[2] = 0.0;
		bcs.push_back(bc);
	}
	else
	{
		int i = n.size() - 6;

		//first param: 1st spline, 2nd spline, 3rd spline
		bc[0] = -1.0;
		bc[1] = 0.0;
		bc[2] = 1.0;
		bcs.push_back(bc);

		//second param: 1st spline, 2nd spline, 3rd spline
		bc[0] = 2.0;
		bc[1] = 1.0;
		bc[2] = 0.0;
		bcs.push_back(bc);
	}
}

void SetBoundaryConditions3_1D_OR_4(vector<double> n, vector<vector<double> >& bcs, bool uniform)
{
	//first derivatie = 0, second derivative = 0
	vector<double> bc(3);

	if (uniform)
	{
		bc[0] = 1.0;
		bc[1] = 1.0;
		bc[2] = 1.0;
		bcs.push_back(bc);
	}
	else
	{
		int i = n.size() - 6;

		//first param: 1st spline, 2nd spline, 3rd spline
		bc[0] = 1.0;
		bc[1] = 1.0;
		bc[2] = 1.0;
		bcs.push_back(bc);
	}
}

void SetBoundaryConditions3_1D_CO_1(vector<double> n, vector<vector<double> >& bcs, double rc, bool uniform)
{
	//first derivatie = 0
	vector<double> bc(3);

	if (uniform)
	{
		bc[0] = 1.0;
		bc[1] = 0.0;
		bc[2] = 0.0;
		bcs.push_back(bc);

		bc[0] = 0.0;
		bc[1] = 1.0;
		bc[2] = 0.0;
		bcs.push_back(bc);

		bc[0] = 0.0;
		bc[1] = 0.0;
		bc[2] = 1.0;
		bcs.push_back(bc);
	}
	else
	{
		int i = n.size() - 6;

		//-3rd spline: -3rd spline, -2nd spline, last spline
		bc[0] = 1.0;
		bc[1] = 0.0;
		bc[2] = 0.0;
		bcs.push_back(bc);

		//-2nd spline: -3rd spline, -2nd spline, last spline
		bc[0] = 0.0;
		bc[1] = 1.0;
		bc[2] = 0.0;
		bcs.push_back(bc);

		//-2nd spline: -3rd spline, -2nd spline, last spline
		bc[0] = 0.0;
		bc[1] = 0.0;
		bc[2] = 1.0;
		bcs.push_back(bc);
	}
}

void SetBoundaryConditions3_1D_CO_2(vector<double> n, vector<vector<double> >& bcs, double rc, bool uniform)
{
	//first derivatie = 0
	vector<double> bc(3);

	if (uniform)
	{
		bc[0] = 1.0;
		bc[1] = 0.0;
		bc[2] = 1.0;
		bcs.push_back(bc);

		bc[0] = 0.0;
		bc[1] = 1.0;
		bc[2] = 0.0;
		bcs.push_back(bc);
	}
	else
	{
		int i = n.size() - 6;

		//-3rd spline: -3rd spline, -2nd spline, last spline
		bc[0] = 1.0;
		bc[1] = 0.0;
		bc[2] = ((n[i + 3] - n[i + 4]) * (n[i + 2] - n[i + 5])) / ((n[i + 2] - n[i + 3]) * (n[i + 1] - n[i + 4]));
		bcs.push_back(bc);

		//-2nd spline: -3rd spline, -2nd spline, last spline
		bc[0] = 0.0;
		bc[1] = 1.0;
		bc[2] = (n[i + 1] * (n[i + 2] - n[i + 3]) - n[i + 2] * n[i + 3] + n[i + 3] * n[i + 4] + n[i + 3] * n[i + 5] - n[i + 4] * n[i + 5]) / ((n[i + 2] - n[i + 3]) * (n[i + 1] - n[i + 4]));
		bcs.push_back(bc);
	}
}

void SetBoundaryConditions3_1D_CO_3(vector<double> n, vector<vector<double> >& bcs, double rc, bool uniform)
{
	//first derivatie = 0
	vector<double> bc(3);

	if (uniform)
	{
		bc[0] = 1.0;
		bc[1] = 0.0;
		bc[2] = -1.0;
		bcs.push_back(bc);

		bc[0] = 0.0;
		bc[1] = 1.0;
		bc[2] = 2.0;
		bcs.push_back(bc);
	}
	else
	{
		int i = n.size() - 6;

		//-3rd spline: -3rd spline, -2nd spline, last spline
		bc[0] = 1.0;
		bc[1] = 0.0;
		bc[2] = -1.0;
		bcs.push_back(bc);

		//-2nd spline: -3rd spline, -2nd spline, last spline
		bc[0] = 0.0;
		bc[1] = 1.0;
		bc[2] = 2.0;
		bcs.push_back(bc);
	}
}

void SetBoundaryConditions3_1D_CO_4(vector<double> n, vector<vector<double> >& bcs, double rc, bool uniform)
{
	//first derivatie = 0, second derivative = 0
	vector<double> bc(3);

	if (uniform)
	{
		bc[0] = 1.0;
		bc[1] = 1.0;
		bc[2] = 1.0;
		bcs.push_back(bc);
	}
	else
	{
		int i = n.size() - 6;

		//-3rd spline: -3rd spline, -2nd spline, last spline
		bc[0] = 1.0;
		bc[1] = 1.0;
		bc[2] = 1.0;
		bcs.push_back(bc);
	}
}

void SetBoundaryConditions4_1D_OR_1(vector<double> n, vector<vector<double> >& bcs)
{
	vector<double> bc(4);

	int i = n.size() - 8;

	/*
	 bc[0] = 1.0;
	 bc[1] = 1.0;
	 bc[2] = 0.0;
	 bc[3] = 0.0;
	 bcs.push_back(bc);

	 bc[0] = 1.0;
	 bc[1] = 0.0;
	 bc[2] = 1.0;
	 bc[3] = 0.0;
	 bcs.push_back(bc);

	 bc[0] = -1.0;
	 bc[1] = 0.0;
	 bc[2] = 0.0;
	 bc[3] = 1.0;
	 bcs.push_back(bc);
	 */

	bc[0] = 1.5;
	bc[1] = 0.5;
	bc[2] = 1.0;
	bc[3] = 0.0;
	bcs.push_back(bc);

	bc[0] = -0.5;
	bc[1] = 0.5;
	bc[2] = 0.0;
	bc[3] = 1.0;
	bcs.push_back(bc);
}

void SetBoundaryConditions4_1D_CO_2(vector<double> n, vector<vector<double> >& bcs, double rc)
{
	vector<double> bc(4);

	int i = n.size() - 8;

	bc[0] = 1.0;
	bc[1] = 0.0;
	bc[2] = 0.5;
	bc[3] = -0.5;
	bcs.push_back(bc);

	bc[0] = 0.0;
	bc[1] = 1.0;
	bc[2] = 0.5;
	bc[3] = 1.5;
	bcs.push_back(bc);
}

}

