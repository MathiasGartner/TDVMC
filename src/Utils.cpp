#include "Utils.h"

using namespace std;

bool FileExist(string path)
{
	ifstream file(path.c_str());
	return file.good();
}

void RemoveFile(string path)
{
	remove(path.c_str());
}

bool CopyFile(string filePathSource, string filePathNew)
{
	ifstream fileSource(filePathSource, ifstream::binary);
	ofstream fileNew(filePathNew, ifstream::binary);
	fileNew << fileSource.rdbuf();
	return fileNew && fileSource;
}

string PrintArrayValues(double* r, int length)
{
	string s = "";
	for (int i = 0; i < length; i++)
	{
		s = s + to_string(r[i]) + ", ";
	}
	return s;
}

double VectorDotProduct(vector<double>& v1, vector<double>& v2)
{
	double sum = 0;
	for (unsigned int i = 0; i < v1.size(); i++)
	{
		sum += v1[i] * v2[i];
	}
	return sum;
}

double VectorNorm2(vector<double>& r)
{
	double sum = 0;
	for (unsigned int i = 0; i < r.size(); i++)
	{
		sum += pow(r[i], 2);
	}
	return sum;
}

double VectorNorm(vector<double>& r)
{
	double norm;
	norm = sqrt(VectorNorm2(r));
	return norm;
}

double VectorDisplacement(vector<double>& ri, vector<double>& rj, vector<double>& result)
{
	double norm;
	result.resize(ri.size());
	for (unsigned int i = 0; i < ri.size(); i++)
	{
		result[i] = ri[i] - rj[i];
	}
	norm = VectorNorm(result);
	return norm;
}

double GetCoordinateNIC(double r)
{
	return r - LBOX * round(r / LBOX);
}

void GetVectorNIC(vector<double>& r, vector<double>& rNIC)
{
	rNIC.resize(r.size());
	for (unsigned int i = 0; i < r.size(); i++)
	{
		rNIC[i] = GetCoordinateNIC(r[i]);
	}
}

void VectorDiffNIC(vector<double>& ri, vector<double>& rj, vector<double>& result)
{
	double delta;
	result.resize(ri.size());
	for (unsigned int i = 0; i < ri.size(); i++)
	{
		delta = ri[i] - rj[i];
		result[i] = GetCoordinateNIC(delta);
	}
}

double VectorDisplacementNIC(vector<double>& ri, vector<double>& rj, vector<double>& result)
{
	double norm;
	VectorDiffNIC(ri, rj, result);
	norm = VectorNorm(result);
	return norm;
}

double Mean(vector<double>& v)
{
	double mean;
	int nvalues = v.size();
	double sum = 0;
	for (int i = 0; i < nvalues; i++)
	{
		sum += v[i];
	}
	mean = sum / (double) nvalues;
	return mean;
}

vector<double> Mean(vector<vector<double> >& v)
{
	int length = v.size();
	int nvalues = v[0].size();
	vector<double> means(nvalues);
	vector<double> sums(nvalues, 0);
	for (int i = 0; i < length; i++)
	{
		for (int j = 0; j < nvalues; j++)
		{
			sums[j] += v[i][j];
		}
	}
	for (int i = 0; i < nvalues; i++)
	{
		means[i] = sums[i] / length;
	}
	return means;
}

vector<vector<double> > Mean(vector<vector<vector<double> > >& v)
{
	int length = v.size();
	int length2 = v[0].size();
	int nvalues = v[0][0].size();
	vector<vector<double> > means(length2);
	vector<vector<double> > sums(length2);
	for (int i = 0; i < length2; i++)
	{
		means[i].clear();
		means[i].resize(nvalues, 0);
		sums[i].clear();
		sums[i].resize(nvalues, 0);
	}
	for (int i = 0; i < length; i++)
	{
		for (int j = 0; j < length2; j++)
		{
			for (int m = 0; m < nvalues; m++)
			{
				sums[j][m] += v[i][j][m];
			}
		}
	}
	for (int i = 0; i < length2; i++)
	{
		for (int j = 0; j < nvalues; j++)
		{
			means[i][j] = sums[i][j] / (double) length;
		}
	}
	return means;
}

//vector<double> Mean(vector<vector<double> >* v)
//{
//	vector<double> mean;
//	int nvalues = (*v).size();
//	for (int i = 0; i < nvalues; i++)
//	{
//		mean.push_back(Mean(&((*v)[i])));
//	}
//	return mean;
//}
//
//vector<vector<double> > Mean(vector<vector<vector<double> > >* v)
//{
//	vector<vector<double> > mean;
//	int nvalues = (*v).size();
//	for (int i = 0; i < nvalues; i++)
//	{
//		mean.push_back(Mean(&((*v)[i])));
//	}
//	return mean;
//}

void CopyArray(double* original, double* copy, int dim)
{
	for (int i = 0; i < dim; i++)
	{
		copy[i] = original[i];
	}
}

void Copy2DArray(double** original, double** copy, int dim1, int dim2)
{
	for (int i = 0; i < dim1; i++)
	{
		for (int j = 0; j < dim2; j++)
		{
			copy[i][j] = original[i][j];
		}
	}
}

vector<double> ArrayToVector(double* arr, int length)
{
	vector<double> vec;
	for (int i = 0; i < length; i++)
	{
		vec.push_back(arr[i]);
	}
	return vec;
}

void split(const string &s, char delim, vector<string> &elems)
{
	stringstream ss;
	ss.str(s);
	string item;
	while (getline(ss, item, delim))
	{
		elems.push_back(item);
	}
}

vector<string> split(const string &s, char delim)
{
	vector<string> elems;
	split(s, delim, elems);
	return elems;
}


double GetSizeOfVector(vector<double>& v)
{
	return (sizeof(v) + sizeof(double) * v.capacity()) / 1024.0 / 1024.0;
}

double GetSizeOfVector(vector<vector<double> >& v)
{
	double size = 0;
	for (auto &element : v)
	{
		size += GetSizeOfVector(element);
	}
	return size;
}

double GetSizeOfVector(vector<vector<vector<double> > >& v)
{
	double size = 0;
	for (auto &element : v)
	{
		size += GetSizeOfVector(element);
	}
	return size;
}

void ClearVector(vector<double>& v)
{
	fill(v.begin(), v.end(), 0);
}

void ClearVector(vector<vector<double> >& v)
{
	for (auto &i : v)
	{
		ClearVector(i);
	}
}

void ClearVector(vector<vector<vector<double> > >& v)
{
	for (auto &i : v)
	{
		ClearVector(i);
	}
}

void ClearVector(vector<vector<vector<vector<double> > > >& v)
{
	for (auto &i : v)
	{
		ClearVector(i);
	}
}

void WriteDataToFile(double* data, int n, string filename, string header)
{
	ofstream file;
	//cout << OUTPUT_DIRECTORY << filename << ".csv" << endl;
	file.open(OUTPUT_DIRECTORY + filename + ".csv", ios::out);
	file.precision(exportPrecision);
	file << header << endl;
	for (int i = 0; i < n; i++)
	{
		file << fixed << data[i] << endl;
	}
	file.close();
}

void WriteDataToFile(double** data, int n1, int n2, string filename, string header)
{
	ofstream file;
	//cout << OUTPUT_DIRECTORY << filename << ".csv" << endl;
	file.open(OUTPUT_DIRECTORY + filename + ".csv", ios::out);
	file.precision(exportPrecision);
	file << header << endl;
	for (int i = 0; i < n1; i++)
	{
		for (int j = 0; j < n2; j++)
		{
			file << fixed << data[i][j];
			if (j < n2 - 1)
			{
				file << ", ";
			}
		}
		if (i < n1 - 1)
		{
			file << endl;
		}
	}
	file.close();
}

void WriteDataToFile(double data, string filename, string header)
{
	ofstream file;
	//cout << OUTPUT_DIRECTORY << filename << ".csv" << endl;
	file.open(OUTPUT_DIRECTORY + filename + ".csv", ios::out);
	file.precision(exportPrecision);
	file << header << endl;
	file << fixed << data << endl;
	file.close();
}

void WriteDataToFile(vector<double>& data, string filename, string header, int everyNth)
{
	ofstream file;
	//cout << OUTPUT_DIRECTORY << filename << ".csv" << endl;
	file.open(OUTPUT_DIRECTORY + filename + ".csv", ios::out);
	file.precision(exportPrecision);
	file << header << endl;
	for (unsigned int i = 0; i < data.size(); i += everyNth)
	{
		file << fixed << data[i] << endl;
	}
	file.close();
}

void WriteDataToFile(vector<vector<double> >& data, string filename, string header, int everyNth, bool writeHeader)
{
	ofstream file;
	//cout << OUTPUT_DIRECTORY << filename << ".csv" << endl;
	file.open(OUTPUT_DIRECTORY + filename + ".csv", ios::out);
	file.precision(exportPrecision);
	if (writeHeader)
	{
		file << header << endl;
	}
	for (unsigned int i = 0; i < data.size(); i += everyNth)
	{
		file << fixed << JoinVector(data[i]) << endl;
	}
	file.close();
}

void WriteDataToFile(vector<vector<vector<double> > >& data, string filename, string header)
{
	ofstream file;
	//cout << OUTPUT_DIRECTORY << filename << ".csv" << endl;
	file.open(OUTPUT_DIRECTORY + filename + ".m", ios::out);
	file.precision(exportPrecision);
	//file << header << endl;
	file << "{" << endl;
	for (unsigned int i = 0; i < data.size(); i++)
	{
		file << "{" << endl;
		for (unsigned int j = 0; j < data[i].size(); j++)
		{
			file << "{";
			file << fixed << JoinVector(data[i][j]);
			file << (j < data[i].size() - 1 ? "}," : "}") << endl;
		}
		file << (i < data.size() - 1 ? "}," : "}") << endl;
	}
	file << "}" << endl;
	file.close();
}

vector<vector<vector<double> > > ReadFromFile(string filePath, int headerlines)
{
	vector<vector<vector<double> > > data;

	Json::Value configData;
	Json::Reader configReader;
	ifstream configFile(filePath, ifstream::binary);
	configReader.parse(configFile, configData, false);

	auto allKValues = configData["data"];
	for (unsigned int i = 0; i < allKValues.size(); i++)
	{
		data.push_back(vector<vector<double> >());
		auto allEqualNormVectors = allKValues[i];
		for (unsigned int j = 0; j < allEqualNormVectors.size(); j++)
		{
			data[i].push_back(vector<double>());
			auto kValues = allEqualNormVectors[j];
			for (unsigned int l = 0; l < kValues.size(); l++)
			{
				data[i][j].push_back(kValues[l].asDouble());
			}
		}
	}
	return data;
}

int get_cpu_id()
{
	//INFO: Get the the current process' stat file from the proc filesystem
	int cpu_id = -1;
	FILE* procfile = fopen("/proc/self/stat", "r");
	long to_read = 8192;
	char buffer[to_read];
	int read = fread(buffer, sizeof(char), to_read, procfile);
	if (read > 0)
	{
		fclose(procfile);

		// Field with index 38 (zero-based counting) is the one we want
		char* line = strtok(buffer, " ");
		for (int i = 1; i < 38; i++)
		{
			line = strtok(NULL, " ");
		}

		line = strtok(NULL, " ");
		cpu_id = atoi(line);
	}
	return cpu_id;
}
