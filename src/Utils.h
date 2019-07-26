#pragma once

#include "Constants.h"
#include "MathOperators.h"

#include "Observables/Observable.h"
#include "Observables/ObservableV.h"
#include "Observables/ObservableVsOnGrid.h"

#include <cmath>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <random>
#include <string>
#include <sstream>

#include "../resources/json/json.h"
#include "../resources/json/json-forwards.h"

using namespace std;

const int exportPrecision = 10;
const int exportPrecisionSc = 8;
const int colWidth = 20;
const string fileExtension = ".dat";

template<typename T> string JoinVector(const vector<T>& v, int everyNth = 1, bool fixedForm = false)
{
	string text = "";

	for (unsigned int i = 0; i < v.size(); i += everyNth)
	{
		ostringstream oss;
		oss.precision(exportPrecision);
		oss << (fixedForm ? std::fixed : std::scientific) << v[i];
		text += oss.str() + ",";
	}
	if (text.size() > 0)
	{
		text = text.substr(0, text.size() - 1);
	}
	return text;
}

bool FileExist(string path);
void RemoveFile(string path);
bool CopyFile(string filePathSource, string filePathNew);

string PrintArrayValues(double* r, int length);

extern double (*VectorDotProduct_DIM)(vector<double>&, vector<double>&);
double VectorDotProduct_1D(vector<double>& v1, vector<double>& v2);
double VectorDotProduct_2D(vector<double>& v1, vector<double>& v2);
double VectorDotProduct_3D(vector<double>& v1, vector<double>& v2);
double VectorDotProduct(vector<double>& v1, vector<double>& v2);

extern double (*VectorNorm2_DIM)(vector<double>& r);
double VectorNorm2_1D(vector<double>& r);
double VectorNorm2_2D(vector<double>& r);
double VectorNorm2_3D(vector<double>& r);
double VectorNorm2(vector<double>& r);

extern double (*VectorNorm_DIM)(vector<double>& r);
double VectorNorm_1D(vector<double>& r);
double VectorNorm_2D(vector<double>& r);
double VectorNorm_3D(vector<double>& r);
double VectorNorm(vector<double>& r);

extern double (*VectorDisplacement_DIM)(vector<double>& ri, vector<double>& rj, vector<double>& result);
double VectorDisplacement_1D(vector<double>& ri, vector<double>& rj, vector<double>& result);
double VectorDisplacement_2D(vector<double>& ri, vector<double>& rj, vector<double>& result);
double VectorDisplacement_3D(vector<double>& ri, vector<double>& rj, vector<double>& result);
double VectorDisplacement(vector<double>& ri, vector<double>& rj, vector<double>& result);

double GetCoordinateNIC(double r);

extern void (*GetVectorNIC_DIM)(vector<double>& r, vector<double>& rNIC);
void GetVectorNIC_1D(vector<double>& r, vector<double>& rNIC);
void GetVectorNIC_2D(vector<double>& r, vector<double>& rNIC);
void GetVectorNIC_3D(vector<double>& r, vector<double>& rNIC);
void GetVectorNIC(vector<double>& r, vector<double>& rNIC);

extern void (*VectorDiffNIC_DIM)(vector<double>& ri, vector<double>& rj, vector<double>& result);
void VectorDiffNIC_1D(vector<double>& ri, vector<double>& rj, vector<double>& result);
void VectorDiffNIC_2D(vector<double>& ri, vector<double>& rj, vector<double>& result);
void VectorDiffNIC_3D(vector<double>& ri, vector<double>& rj, vector<double>& result);
void VectorDiffNIC(vector<double>& ri, vector<double>& rj, vector<double>& result);

extern double (*VectorDisplacementNIC_DIM)(vector<double>& ri, vector<double>& rj, vector<double>& result);
double VectorDisplacementNIC_1D(vector<double>& ri, vector<double>& rj, vector<double>& result);
double VectorDisplacementNIC_2D(vector<double>& ri, vector<double>& rj, vector<double>& result);
double VectorDisplacementNIC_3D(vector<double>& ri, vector<double>& rj, vector<double>& result);
double VectorDisplacementNIC(vector<double>& ri, vector<double>& rj, vector<double>& result);

double GetCornerAngle(vector<double>& r1, vector<double>& r2, vector<double>& r3);

vector<double> AllVectorDisplacements(vector<double>& ri, vector<double>& rj, double maxDistance);
vector<double> AllVectorDisplacements(vector<double>& ri, vector<double>& rj, double maxDistance, vector<vector<double> >& eVecs);

double GetRelativeError(double d1, double d2);
double GetRelativeErrorLastElements(vector<double> v);

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

extern void (*AssignVector_DIM)(vector<double>& original, vector<double>& copy);
void AssignVector_1D(vector<double>& original, vector<double>& copy);
void AssignVector_2D(vector<double>& original, vector<double>& copy);
void AssignVector_3D(vector<double>& original, vector<double>& copy);
void AssignVector(vector<double>& original, vector<double>& copy);

double GetSizeOfVector(vector<double>& v);
double GetSizeOfVector(vector<vector<double> >& v);
double GetSizeOfVector(vector<vector<vector<double> > >& v);

void ClearVector(vector<long long>& v);
void ClearVector(vector<double>& v);
void ClearVector(vector<vector<double> >& v);
void ClearVector(vector<vector<vector<double> > >& v);
void ClearVector(vector<vector<vector<vector<double> > > >& v);

void InitVector(vector<int>& v, int size, int initialValue);
void InitVector(vector<vector<int> >& v, int size1, int size2, int initialValue);
void InitVector(vector<long long>& v, int size, long long initialValue);
void InitVector(vector<double>& v, int size, double initialValue);
void InitVector(vector<vector<double> >& v, int size1, int size2, double initialValue);
void InitVector(vector<vector<vector<double> > >& v, int size1, int size2, int size3, double initialValue);

void SetFileFormat(ofstream& f);
void WriteLineToFile(ofstream& f, double data);
void WriteLineToFile(ofstream& f, vector<double>& data);

void WriteDataToFile(double data, string filename, string header);
void WriteDataToFile(vector<double>& data, string filename, string header, int everyNth = 1);
void WriteDataToFile(vector<vector<double> >& data, string filename, string header, int everyNth = 1, bool writeHeader = true);
void WriteDataToFile(vector<vector<vector<double> > >& data, string filename, string header);

void WriteDataToFile(Observables::ObservableCollection& data, string filename);

void AppendDataToFile(double data, string filename);
void AppendDataToFile(vector<double>& data, string filename);

vector<vector<vector<double> > > ReadKValuesFromJsonFile(string filePath);

int get_cpu_id();
