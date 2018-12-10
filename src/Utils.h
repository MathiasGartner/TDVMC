#pragma once

#include "Constants.h"

#include <cmath>
#include <cstring>
#include <fstream>
#include <random>
#include <string>
#include <sstream>

#include "../resources/json/json.h"
#include "../resources/json/json-forwards.h"

using namespace std;

const int exportPrecision = 25;

template<typename T> string JoinVector(const vector<T>& v, int everyNth = 1)
{
	string text = "";

	for (unsigned int i = 0; i < v.size(); i += everyNth)
	{
		ostringstream oss;
		oss.precision(exportPrecision);
		oss << fixed << v[i];
		text += oss.str() + ",";
	}
	if (text.size() > 0)
	{
		text = text.substr(0, text.size() - 2);
	}
	return text;
}

bool FileExist(string path);
void RemoveFile(string path);
bool CopyFile(string filePathSource, string filePathNew);

string PrintArrayValues(double* r, int length);

double VectorDotProduct(vector<double>& v1, vector<double>& v2);
double VectorNorm2(vector<double>& r);
double VectorNorm(vector<double>& r);
double VectorDisplacement(vector<double>& ri, vector<double>& rj, vector<double>& result);

double GetCoordinateNIC(double r);
void GetVectorNIC(vector<double>& r, vector<double>& rNIC);
void VectorDiffNIC(vector<double>& ri, vector<double>& rj, vector<double>& result);
double VectorDisplacementNIC(vector<double>& ri, vector<double>& rj, vector<double>& result);

double Sum(vector<double>& v);
vector<double> Sum(vector<vector<double> >& v);
vector<double> OuterSum(vector<vector<double> >& v);
double Mean(vector<double>& v);
vector<double> Mean(vector<vector<double> >& v);
vector<vector<double> > Mean(vector<vector<vector<double> > >& v);

void CopyArray(double* original, double* copy, int dim);
void Copy2DArray(double** original, double** copy, int dim1, int dim2);
vector<double> ArrayToVector(double* arr, int length);

void split(const string &s, char delim, vector<string> &elems);
vector<string> split(const string &s, char delim);

double GetSizeOfVector(vector<double>& v);
double GetSizeOfVector(vector<vector<double> >& v);
double GetSizeOfVector(vector<vector<vector<double> > >& v);

void ClearVector(vector<double>& v);
void ClearVector(vector<vector<double> >& v);
void ClearVector(vector<vector<vector<double> > >& v);
void ClearVector(vector<vector<vector<vector<double> > > >& v);

void WriteDataToFile(double* data, int n, string filename, string header);
void WriteDataToFile(double** data, int n1, int n2, string filename, string header);
void WriteDataToFile(double data, string filename, string header);
void WriteDataToFile(vector<double>& data, string filename, string header, int everyNth = 1);
void WriteDataToFile(vector<vector<double> >& data, string filename, string header, int everyNth = 1, bool writeHeader = true);
void WriteDataToFile(vector<vector<vector<double> > >& data, string filename, string header);

void AppendDataToFile(vector<double>& data, string filename);

vector<vector<vector<double> > > ReadKValuesFromJsonFile(string filePath);

int get_cpu_id();
