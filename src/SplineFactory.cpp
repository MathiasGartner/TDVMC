
#include "SplineFactory.h"

using namespace std;

namespace SplineFactory
{

vector<vector<vector<double> > > GetWeights(vector<double> n)
{
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
		weights[s][0][0] = pow(n[i],3)/((n[i] - n[i+1])*(n[i] - n[i+2])*(n[i] - n[i+3]));
		weights[s][0][1] = (-3*pow(n[i],2))/((n[i] - n[i+1])*(n[i] - n[i+2])*(n[i] - n[i+3]));
		weights[s][0][2] = (3*n[i])/((n[i] - n[i+1])*(n[i] - n[i+2])*(n[i] - n[i+3]));
		weights[s][0][3] = -(1/((n[i] - n[i+1])*(n[i] - n[i+2])*(n[i] - n[i+3])));
		//Part 2:
		weights[s][1][0] = (pow(n[i],2)*pow(n[i+1],2)*n[i+2] + pow(n[i],2)*pow(n[i+1],2)*n[i+3] + pow(n[i],2)*pow(n[i+1],2)*n[i+4] - pow(n[i],2)*n[i+1]*n[i+3]*n[i+4] - n[i]*pow(n[i+1],2)*n[i+3]*n[i+4] + pow(n[i],2)*n[i+2]*n[i+3]*n[i+4] + n[i]*n[i+1]*n[i+2]*n[i+3]*n[i+4] + pow(n[i+1],2)*n[i+2]*n[i+3]*n[i+4] - pow(n[i],2)*n[i+1]*n[i+2]*(n[i+3] + n[i+4]) - n[i]*pow(n[i+1],2)*n[i+2]*(n[i+3] + n[i+4]))/((n[i+1] - n[i+2])*(-n[i] + n[i+2])*(n[i] - n[i+3])*(n[i+1] - n[i+3])*(n[i+1] - n[i+4]));
		weights[s][1][1] = (-3*pow(n[i],2)*pow(n[i+1],2) + 3*n[i]*n[i+1]*n[i+2]*n[i+3] + 3*n[i]*n[i+1]*n[i+2]*n[i+4] + 3*n[i]*n[i+1]*n[i+3]*n[i+4] - 3*n[i]*n[i+2]*n[i+3]*n[i+4] - 3*n[i+1]*n[i+2]*n[i+3]*n[i+4])/((n[i+1] - n[i+2])*(-n[i] + n[i+2])*(n[i] - n[i+3])*(n[i+1] - n[i+3])*(n[i+1] - n[i+4]));
		weights[s][1][2] = (3*pow(n[i],2)*n[i+1] + 3*n[i]*pow(n[i+1],2) - 3*n[i]*n[i+1]*n[i+2] - 3*n[i]*n[i+1]*n[i+3] - 3*n[i]*n[i+1]*n[i+4] + 3*n[i+2]*n[i+3]*n[i+4])/((n[i+1] - n[i+2])*(-n[i] + n[i+2])*(n[i] - n[i+3])*(n[i+1] - n[i+3])*(n[i+1] - n[i+4]));
		weights[s][1][3] = (-pow(n[i],2) - n[i]*n[i+1] - pow(n[i+1],2) + n[i]*n[i+2] + n[i+1]*n[i+2] - n[i+2]*n[i+3] - n[i+2]*n[i+4] - n[i+3]*n[i+4] + n[i]*(n[i+3] + n[i+4]) + n[i+1]*(n[i+3] + n[i+4]))/((n[i+1] - n[i+2])*(-n[i] + n[i+2])*(n[i] - n[i+3])*(n[i+1] - n[i+3])*(n[i+1] - n[i+4]));
		//Part 3:
		weights[s][2][0] = (n[i]*n[i+1]*n[i+2]*pow(n[i+3],2) + n[i]*n[i+1]*n[i+2]*n[i+3]*n[i+4] - n[i]*n[i+1]*pow(n[i+3],2)*n[i+4] - n[i]*n[i+2]*pow(n[i+3],2)*n[i+4] - n[i+1]*n[i+2]*pow(n[i+3],2)*n[i+4] + n[i]*n[i+1]*n[i+2]*pow(n[i+4],2) - n[i]*n[i+1]*n[i+3]*pow(n[i+4],2) - n[i]*n[i+2]*n[i+3]*pow(n[i+4],2) - n[i+1]*n[i+2]*n[i+3]*pow(n[i+4],2) + n[i]*pow(n[i+3],2)*pow(n[i+4],2) + n[i+1]*pow(n[i+3],2)*pow(n[i+4],2) + n[i+2]*pow(n[i+3],2)*pow(n[i+4],2))/((n[i+2] - n[i+3])*(-n[i] + n[i+3])*(-n[i+1] + n[i+3])*(n[i+1] - n[i+4])*(n[i+2] - n[i+4]));
		weights[s][2][1] = (-3*n[i]*n[i+1]*n[i+2]*n[i+3] - 3*n[i]*n[i+1]*n[i+2]*n[i+4] + 3*n[i]*n[i+1]*n[i+3]*n[i+4] + 3*n[i]*n[i+2]*n[i+3]*n[i+4] + 3*n[i+1]*n[i+2]*n[i+3]*n[i+4] - 3*pow(n[i+3],2)*pow(n[i+4],2))/((n[i+2] - n[i+3])*(-n[i] + n[i+3])*(-n[i+1] + n[i+3])*(n[i+1] - n[i+4])*(n[i+2] - n[i+4]));
		weights[s][2][2] = (3*n[i]*n[i+1]*n[i+2] - 3*n[i]*n[i+3]*n[i+4] - 3*n[i+1]*n[i+3]*n[i+4] - 3*n[i+2]*n[i+3]*n[i+4] + 3*pow(n[i+3],2)*n[i+4] + 3*n[i+3]*pow(n[i+4],2))/((n[i+2] - n[i+3])*(-n[i] + n[i+3])*(-n[i+1] + n[i+3])*(n[i+1] - n[i+4])*(n[i+2] - n[i+4]));
		weights[s][2][3] = (-(n[i]*n[i+1]) - n[i]*n[i+2] - n[i+1]*n[i+2] + n[i]*n[i+3] + n[i+1]*n[i+3] + n[i+2]*n[i+3] - pow(n[i+3],2) + n[i]*n[i+4] + n[i+1]*n[i+4] + n[i+2]*n[i+4] - n[i+3]*n[i+4] - pow(n[i+4],2))/((n[i+2] - n[i+3])*(-n[i] + n[i+3])*(-n[i+1] + n[i+3])*(n[i+1] - n[i+4])*(n[i+2] - n[i+4]));
		//Part 4:
		weights[s][3][0] = -(pow (n[i+4],3)/((n[i+3] - n[i+4])*(-n[i+1] + n[i+4])*(-n[i+2] + n[i+4])));
		weights[s][3][1] = (3*pow(n[i+4],2))/((n[i+3] - n[i+4])*(-n[i+1] + n[i+4])*(-n[i+2] + n[i+4]));
		weights[s][3][2] = (-3*n[i+4])/((n[i+3] - n[i+4])*(-n[i+1] + n[i+4])*(-n[i+2] + n[i+4]));
		weights[s][3][3] = 1/((n[i+3] - n[i+4])*(-n[i+1] + n[i+4])*(-n[i+2] + n[i+4]));
	}
	return weights;
}

vector<double> GetBoundaryConditions1_1(vector<double> n)
{
	vector<double> bc(4);

	int i = n.size() - 6;

	bc[0] = 1.0;
	bc[1] = ((n[i+2] - n[i+3])*(-n[i+1] + n[i+3])*(n[i+3] - n[i+4])*(n[i+2] - n[i+5])*(n[i+2] + n[i+3] - n[i+4] - n[i+5]))/((n[i+1] - n[i+4])*(n[i+2] - n[i+4])*(-((n[i+2] - n[i+3])*n[i+4]*(n[i+2] + n[i+4])) + n[i+4]*(2*n[i+2] - 3*n[i+3] + n[i+4])*n[i+5] + (n[i+3] - n[i+4])*pow(n[i+5],2) + n[i+1]*(pow(n[i+2],2) - pow(n[i+4],2) + n[i+2]*(-n[i+3] + n[i+4] - 2*n[i+5]) + (n[i+3] + n[i+4])*n[i+5])));
	bc[2] = ((n[i+2] - n[i+3])*(n[i+3] - n[i+4])*((n[i+2] - n[i+3])*(n[i+3] - 2*n[i+4])*n[i+4] + (n[i+3]*(n[i+2] + n[i+3]) - 4*n[i+3]*n[i+4] + 2*pow(n[i+4],2))*n[i+5] + n[i+1]*((n[i+3] - 2*n[i+4])*n[i+4] + n[i+3]*n[i+5] - n[i+2]*(2*n[i+3] - 3*n[i+4] + n[i+5]))))/((n[i+1] - n[i+4])*(n[i+2] - n[i+5])*(-((n[i+2] - n[i+3])*n[i+4]*(n[i+2] + n[i+4])) + n[i+4]*(2*n[i+2] - 3*n[i+3] + n[i+4])*n[i+5] + (n[i+3] - n[i+4])*pow(n[i+5],2) + n[i+1]*(pow(n[i+2],2) - pow(n[i+4],2) + n[i+2]*(-n[i+3] + n[i+4] - 2*n[i+5]) + (n[i+3] + n[i+4])*n[i+5])));
	bc[3] = 0.0;

	return bc;
}

}

