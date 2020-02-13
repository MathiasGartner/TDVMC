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

vector<double> VectorAbs(vector<double> v)
{
	for (unsigned int i = 0; i < v.size(); i++)
	{
		v[i] = abs(v[i]);
	}
	return v;
}

double (*VectorDotProduct_DIM)(vector<double>&, vector<double>&);

double VectorDotProduct_1D(vector<double>& v1, vector<double>& v2)
{
	double sum = 0;
	sum = v1[0] * v2[0];
	return sum;
}

double VectorDotProduct_2D(vector<double>& v1, vector<double>& v2)
{
	double sum = 0;
	sum = v1[0] * v2[0] + v1[1] * v2[1];
	return sum;
}

double VectorDotProduct_3D(vector<double>& v1, vector<double>& v2)
{
	double sum = 0;
	sum = v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2];
	return sum;
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

double (*VectorNorm2_DIM)(vector<double>& r);

double VectorNorm2_1D(vector<double>& r)
{
	double sum = 0;
	sum = r[0] * r[0];
	return sum;
}

double VectorNorm2_2D(vector<double>& r)
{
	double sum = 0;
	sum = r[0] * r[0] + r[1] * r[1];
	return sum;
}

double VectorNorm2_3D(vector<double>& r)
{
	double sum = 0;
	sum = r[0] * r[0] + r[1] * r[1] + r[2] * r[2];
	return sum;
}

double VectorNorm2(vector<double>& r)
{
	double sum = 0;
	for (unsigned int i = 0; i < r.size(); i++)
	{
		sum += r[i] * r[i];
	}
	return sum;
}

double (*VectorNorm_DIM)(vector<double>& r);

double VectorNorm_1D(vector<double>& r)
{
	double norm;
	norm = sqrt(VectorNorm2_1D(r));
	return norm;
}

double VectorNorm_2D(vector<double>& r)
{
	double norm;
	norm = sqrt(VectorNorm2_2D(r));
	return norm;
}

double VectorNorm_3D(vector<double>& r)
{
	double norm;
	norm = sqrt(VectorNorm2_3D(r));
	return norm;
}

double VectorNorm(vector<double>& r)
{
	double norm;
	norm = sqrt(VectorNorm2(r));
	return norm;
}

double (*VectorDisplacement_DIM)(vector<double>& ri, vector<double>& rj, vector<double>& result);

double VectorDisplacement_1D(vector<double>& ri, vector<double>& rj, vector<double>& result)
{
	//INFO: result needs to have correct size
	double norm;
	result[0] = ri[0] - rj[0];
	norm = VectorNorm_1D(result);
	return norm;
}

double VectorDisplacement_2D(vector<double>& ri, vector<double>& rj, vector<double>& result)
{
	//INFO: result needs to have correct size
	double norm;
	result[0] = ri[0] - rj[0];
	result[1] = ri[1] - rj[1];
	norm = VectorNorm_2D(result);
	return norm;
}

double VectorDisplacement_3D(vector<double>& ri, vector<double>& rj, vector<double>& result)
{
	//INFO: result needs to have correct size
	double norm;
	result[0] = ri[0] - rj[0];
	result[1] = ri[1] - rj[1];
	result[2] = ri[2] - rj[2];
	norm = VectorNorm_3D(result);
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
	int k;
	k = (int) (r * LBOX_R + ((r >= 0.0) ? 0.5 : -0.5));
	return r - k * LBOX;
}

void (*GetVectorNIC_DIM)(vector<double>& r, vector<double>& rNIC);

void GetVectorNIC_1D(vector<double>& r, vector<double>& rNIC)
{
	rNIC[0] = GetCoordinateNIC(r[0]);
}

void GetVectorNIC_2D(vector<double>& r, vector<double>& rNIC)
{
	rNIC[0] = GetCoordinateNIC(r[0]);
	rNIC[1] = GetCoordinateNIC(r[1]);
}

void GetVectorNIC_3D(vector<double>& r, vector<double>& rNIC)
{
	rNIC[0] = GetCoordinateNIC(r[0]);
	rNIC[1] = GetCoordinateNIC(r[1]);
	rNIC[2] = GetCoordinateNIC(r[2]);
}

void GetVectorNIC(vector<double>& r, vector<double>& rNIC)
{
	for (unsigned int i = 0; i < r.size(); i++)
	{
		rNIC[i] = GetCoordinateNIC(r[i]);
	}
}

void (*VectorDiffNIC_DIM)(vector<double>& ri, vector<double>& rj, vector<double>& result);

void VectorDiffNIC_1D(vector<double>& ri, vector<double>& rj, vector<double>& result)
{
	double delta;
	delta = ri[0] - rj[0];
	result[0] = GetCoordinateNIC(delta);
}

void VectorDiffNIC_2D(vector<double>& ri, vector<double>& rj, vector<double>& result)
{
	double delta;
	delta = ri[0] - rj[0];
	result[0] = GetCoordinateNIC(delta);
	delta = ri[1] - rj[1];
	result[1] = GetCoordinateNIC(delta);
}

void VectorDiffNIC_3D(vector<double>& ri, vector<double>& rj, vector<double>& result)
{
	double delta;
	delta = ri[0] - rj[0];
	result[0] = GetCoordinateNIC(delta);
	delta = ri[1] - rj[1];
	result[1] = GetCoordinateNIC(delta);
	delta = ri[2] - rj[2];
	result[2] = GetCoordinateNIC(delta);
}

void VectorDiffNIC(vector<double>& ri, vector<double>& rj, vector<double>& result)
{
	double delta;
	for (unsigned int i = 0; i < ri.size(); i++)
	{
		delta = ri[i] - rj[i];
		result[i] = GetCoordinateNIC(delta);
	}
}

double (*VectorDisplacementNIC_DIM)(vector<double>& ri, vector<double>& rj, vector<double>& result);

double VectorDisplacementNIC_1D(vector<double>& ri, vector<double>& rj, vector<double>& result)
{
	double norm;
	VectorDiffNIC_1D(ri, rj, result);
	norm = VectorNorm_1D(result);
	return norm;
}

double VectorDisplacementNIC_2D(vector<double>& ri, vector<double>& rj, vector<double>& result)
{
	double norm;
	VectorDiffNIC_2D(ri, rj, result);
	norm = VectorNorm_2D(result);
	return norm;
}

double VectorDisplacementNIC_3D(vector<double>& ri, vector<double>& rj, vector<double>& result)
{
	double norm;
	VectorDiffNIC_3D(ri, rj, result);
	norm = VectorNorm_3D(result);
	return norm;
}

double VectorDisplacementNIC(vector<double>& ri, vector<double>& rj, vector<double>& result)
{
	double norm;
	VectorDiffNIC(ri, rj, result);
	norm = VectorNorm(result);
	return norm;
}

double GetCornerAngle(vector<double>& r1, vector<double>& r2, vector<double>& r3)
{
	double angle = 0;
	vector<double> tmp(DIM);
	double r12 = VectorDisplacement(r1, r2, tmp);
	double r13 = VectorDisplacement(r1, r3, tmp);
	double r23 = VectorDisplacement(r2, r3, tmp);

	//arccos(a^2 + b^2 - c^2 / 2 a b)
	angle = acos((r12 * r12 + r23 * r23 - r13 * r13) / (2 * r12 * r23));
	angle = angle / M_PI * 180.0; //INFO: in degree
	return angle;
}

vector<double> AllVectorDisplacements(vector<double>& ri, vector<double>& rj, double maxDistance)
{
	//TODO: if i==j the distances are simply L, sqrt(2 L), 2L, ... for 2D
	//TODO: do not include term with i=j=0 in periodics. this term is already calculated by VectorDiffNIC(ri, rj, originalDiff)
	vector<double> displacements;
	vector<double> originalDiff(ri.size());
	VectorDiffNIC(ri, rj, originalDiff);
	vector<vector<int> > periodics;
	int max = 1;
	for (int i = -max; i <= max; i++)
	{
		for (int j = -max; j <= max; j++)
		{
			periodics.push_back(vector<int>( { i, j }));
		}
	}
	double newDistance;
	vector<double> newVec;
	for (unsigned int i = 0; i < periodics.size(); i++)
	{
		newVec = originalDiff + (periodics[i] * LBOX);
		newDistance = VectorNorm(newVec);
		if (newDistance < maxDistance && newDistance > 0)
		{
			displacements.push_back(newDistance);
		}
	}
	return displacements;
}

vector<double> AllVectorDisplacements(vector<double>& ri, vector<double>& rj, double maxDistance, vector<vector<double> >& eVecs)
{
	vector<double> displacements;
	vector<double> originalDiff(ri.size());
	VectorDiffNIC(ri, rj, originalDiff);
	vector<vector<int> > periodics;
	int max = 1;
	for (int i = -max; i <= max; i++)
	{
		for (int j = -max; j <= max; j++)
		{
			periodics.push_back(vector<int>( { i, j }));
		}
	}
	double newDistance;
	vector<double> newVec;
	for (unsigned int i = 0; i < periodics.size(); i++)
	{
		newVec = originalDiff + (periodics[i] * LBOX);
		newDistance = VectorNorm(newVec);
		if (newDistance < maxDistance && newDistance > 0)
		{
			displacements.push_back(newDistance);
			eVecs.push_back(newVec / newDistance);
		}
	}
	return displacements;
}

double GetRelativeError(double d1, double d2)
{
	return 1.0 - d2 / d1;
}

double GetRelativeErrorLastElements(vector<double> v)
{
	return GetRelativeError(v[v.size() - 2], v[v.size() - 1]);
}

double Sum(vector<double>& v)
{
	double sum = 0;
	for (double d : v)
	{
		sum += d;
	}
	return sum;
}

vector<double> Sum(vector<vector<double> >& v)
{
	vector<double> sums;
	for (auto w : v)
	{
		sums.push_back(Sum(w));
	}
	return sums;
}

vector<double> OuterSum(vector<vector<double> >& v)
{
	vector<double> sums;
	double sum;
	for (unsigned int i = 0; i < v[0].size(); i++)
	{
		sum = 0;
		for (unsigned int j = 0; j < v.size(); j++)
		{
			sum += v[j][i];
		}
		sums.push_back(sum);
	}
	return sums;
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

void (*AssignVector_DIM)(vector<double>& original, vector<double>& copy);

void AssignVector_1D(vector<double>& original, vector<double>& copy)
{
	copy[0] = original[0];
}

void AssignVector_2D(vector<double>& original, vector<double>& copy)
{
	copy[0] = original[0];
	copy[1] = original[1];
}

void AssignVector_3D(vector<double>& original, vector<double>& copy)
{
	copy[0] = original[0];
	copy[1] = original[1];
	copy[2] = original[2];
}

void AssignVector(vector<double>& original, vector<double>& copy)
{
	for (unsigned int i = 0; i < original.size(); i++)
	{
		copy[i] = original[i];
	}
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

void ClearVector(vector<long long>& v)
{
	fill(v.begin(), v.end(), 0);
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

//TODO: use templates here!
void InitVector(vector<int>& v, int size, int initialValue)
{
	v.resize(size);
	fill(v.begin(), v.end(), initialValue);
}

void InitVector(vector<vector<int> >& v, int size1, int size2, int initialValue)
{
	v.resize(size1);
	for (auto& i : v)
	{
		InitVector(i, size2, initialValue);
	}
}

void InitVector(vector<long long>& v, int size, long long initialValue)
{
	v.resize(size);
	fill(v.begin(), v.end(), initialValue);
}

void InitVector(vector<double>& v, int size, double initialValue)
{
	v.resize(size);
	fill(v.begin(), v.end(), initialValue);
}

void InitVector(vector<vector<double> >& v, int size1, int size2, double initialValue)
{
	v.resize(size1);
	for (auto& i : v)
	{
		InitVector(i, size2, initialValue);
	}
}

void InitVector(vector<vector<vector<double> > >& v, int size1, int size2, int size3, double initialValue)
{
	v.resize(size1);
	for (auto& i : v)
	{
		InitVector(i, size2, size3, initialValue);
	}
}

void SetFileFormat(ofstream& f)
{
	f.precision(exportPrecisionSc);
	f << std::scientific;
	f << left;
}

void WriteLineToFile(ofstream& f, double data)
{
	f << setw(colWidth) << data;
	f << endl;
}

void WriteLineToFile(ofstream& f, vector<double>& data)
{
	for (auto& value : data)
	{
		f << setw(colWidth) << value;
	}
	f << endl;
}

void WriteDataToFile(double data, string filename, string header)
{
	ofstream file;
	//cout << OUT_DIR << filename << fileExtension << endl;
	file.open(OUT_DIR + filename + fileExtension, ios::out);
	SetFileFormat(file);
	file << header << endl;
	WriteLineToFile(file, data);
	file.close();
}

void WriteDataToFile(vector<double>& data, string filename, string header, int everyNth)
{
	ofstream file;
	//cout << OUT_DIR << filename << fileExtension << endl;
	file.open(OUT_DIR + filename + fileExtension, ios::out);
	SetFileFormat(file);
	file << header << endl;
	for (unsigned int i = everyNth - 1; i < data.size(); i += everyNth)
	{
		WriteLineToFile(file, data[i]);
	}
	file.close();
}

void WriteDataToFile(vector<vector<double> >& data, string filename, string header, int everyNth, bool writeHeader)
{
	ofstream file;
	//cout << OUT_DIR << filename << fileExtension << endl;
	file.open(OUT_DIR + filename + fileExtension, ios::out);
	SetFileFormat(file);
	if (writeHeader)
	{
		file << header << endl;
	}
	for (unsigned int i = everyNth - 1; i < data.size(); i += everyNth)
	{
		WriteLineToFile(file, data[i]);
	}
	file.close();
}

void WriteDataToFile(vector<vector<vector<double> > >& data, string filename, string header)
{
	//TODO: also use scientific notation
	ofstream file;
	//cout << OUT_DIR << filename << ".m" << endl;
	file.open(OUT_DIR + filename + ".m", ios::out);
	file.precision(exportPrecision);
	//file << header << endl;
	file << "{" << endl;
	for (unsigned int i = 0; i < data.size(); i++)
	{
		file << "{" << endl;
		for (unsigned int j = 0; j < data[i].size(); j++)
		{
			file << "{";
			file << std::fixed << JoinVector(data[i][j]);
			file << (j < data[i].size() - 1 ? "}," : "}") << endl;
		}
		file << (i < data.size() - 1 ? "}," : "}") << endl;
	}
	file << "}" << endl;
	file.close();
}

void WriteDataToFile(Observables::ObservableCollection& data, string filename)
{
	for (auto& obs : data.observables)
	{
		ofstream file;
		file.open(OUT_DIR + filename + "_" + obs->name + fileExtension, ios::out);
		SetFileFormat(file);

		if (auto o = dynamic_cast<Observables::ObservableVsOnMultiGrid*>(obs))
		{
			for (auto& j : o->grids)
			{
				file << setw(colWidth) << j.name;
			}
			for (auto& j : o->observablesV)
			{
				file << setw(colWidth) << j.name;
			}
			file << endl;

			for (int i = 0; i < o->totalGridPoints; i++)
			{
				int offset = 1;
				for (auto&j : o->grids)
				{
					int gridIndex = i / offset;
					file << setw(colWidth) << j.grid[gridIndex];
					offset *= j.count;
				}
				for (auto& j : o->observablesV)
				{
					file << setw(colWidth) << j.values[i];
				}
				file << endl;
			}
		}
		else if (auto o = dynamic_cast<Observables::ObservableVsOnGrid*>(obs))
		{
			file << setw(colWidth) << o->grid.name;
			for (auto& j : o->observablesV)
			{
				file << setw(colWidth) << j.name;
			}
			file << endl;
			for (int i = 0; i < o->grid.count; i++)
			{
				file << setw(colWidth) << o->grid.grid[i];
				for (auto& j : o->observablesV)
				{
					file << setw(colWidth) << j.values[i];
				}
				file << endl;
			}
		}
		else if (auto o = dynamic_cast<Observables::ObservableV*>(obs))
		{
			string header = o->name;
			file << header << endl;
			for (unsigned int i = 0; i < o->values.size(); i++)
			{
				file << o->values[i] << endl;
			}
		}
		else if (auto o = dynamic_cast<Observables::Observable*>(obs))
		{
			file << o->value << endl;
		}
		else
		{
			throw std::invalid_argument("ReduceToAverage of IObservable* not implemented for given type.");
		}

		file.close();
	}
}

void AppendDataToFile(double data, string filename)
{
	ofstream file;
	//cout << OUT_DIR << filename << fileExtension << endl;
	file.open(OUT_DIR + filename + fileExtension, ios::app);
	SetFileFormat(file);
	WriteLineToFile(file, data);
	file.close();
}

void AppendDataToFile(vector<double>& data, string filename)
{
	ofstream file;
	//cout << OUT_DIR << filename << fileExtension << endl;
	file.open(OUT_DIR + filename + fileExtension, ios::app);
	SetFileFormat(file);
	WriteLineToFile(file, data);
	file.close();
}

vector<vector<vector<double> > > ReadKValuesFromJsonFile(string filePath)
{
	vector<vector<vector<double> > > data;
	string errors;
	bool successful;
	Json::Value configData;
	Json::CharReaderBuilder builder;
	ifstream configFile(filePath, ifstream::binary);

	successful = Json::parseFromStream(builder, configFile, &configData, &errors);
	if (successful)
	{
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
	}
	return data;
}

void ReadDataFromFile(vector<double>& data, string filePath)
{
	ifstream file(filePath);
	double value;
	while (file >> value)
	{
		data.push_back(value);
	}
	file.close();
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
