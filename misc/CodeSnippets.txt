
void BroadcastConfigItem(ConfigItem* ci)
{
	if (ci->type == STRING)
	{
		MPIMethods::BroadcastValue((string*) (ci->variable), 300);
	}
	else if (ci->type == DOUBLE)
	{
		MPI_Bcast(ci->variable, 1, MPI_DOUBLE, rootRank, MPI_COMM_WORLD);
	}
	else if (ci->type == INT)
	{
		MPI_Bcast(ci->variable, 1, MPI_INT, rootRank, MPI_COMM_WORLD);
	}
	else if (ci->type == ARR_DOUBLE)
	{
		vector<double> data = (*(vector<double>*) (ci->variable));
		MPI_Bcast(data.data(), data.size(), MPI_DOUBLE, rootRank, MPI_COMM_WORLD);
	}
}

void ConfigItem::setValue(Json::Value value)
{
	if (value != 0)
	{
		if (this->type == STRING)
		{
			*(string*) (this->variable) = value.asString();
		}
		else if (this->type == DOUBLE)
		{
			*(double*) (this->variable) = value.asDouble();
		}
		else if (this->type == INT)
		{
			*(int*) (this->variable) = value.asInt();
		}
		else if (this->type == ARR_DOUBLE)
		{
			(*(vector<double>*) (this->variable)).clear();
			for (unsigned int i = 0; i < value.size(); i++)
			{
				(*(vector<double>*) (this->variable)).push_back(value[i].asDouble());
			}
		}
	}
}

string ConfigItem::getJsonString()
{
	string s = "";
	if (this->type == STRING)
	{
		s = "\"" + *(string*) (this->variable) + "\"";
	}
	else if (this->type == DOUBLE)
	{
		s = to_string(*(double*) (this->variable));
	}
	else if (this->type == INT)
	{
		s = to_string(*(int*) (this->variable));
	}
	else if (this->type == ARR_DOUBLE)
	{
		vector<double> v = (*(vector<double>*) (this->variable));
		s = "[ " + JoinVector(v) + " ]";
	}
	return s;
}

double random01()
{
	return ((double) rand() / RAND_MAX);
}

void ReadConfig(string filePath)
{
	//cout << "file exists: " << FileExist(filePath) << "." << endl;
	Json::Value configData;
	Json::Reader configReader;
	ifstream configFile(filePath, ifstream::binary);
	configReader.parse(configFile, configData, false);

	SYSTEM_TYPE = configData["SYSTEM_TYPE"] == 0 ? SYSTEM_TYPE : configData["SYSTEM_TYPE"].asString();
	OUTPUT_DIRECTORY = configData["OUTPUT_DIRECTORY"] == 0 ? OUTPUT_DIRECTORY : configData["OUTPUT_DIRECTORY"].asString();
	N = !configData["N"] ? N : configData["N"].asInt();
	DIM = !configData["DIM"] ? DIM : configData["DIM"].asInt();
	N_PARAM = !configData["N_PARAM"] ? N_PARAM : configData["N_PARAM"].asInt();
	RHO = !configData["RHO"] ? RHO : configData["RHO"].asDouble();
	RC = !configData["RC"] ? RC : configData["RC"].asDouble();
	MC_STEP = !configData["MC_STEP"] ? MC_STEP : configData["MC_STEP"].asDouble();
	MC_STEP_OFFSET = !configData["MC_STEP_OFFSET"] ? MC_STEP_OFFSET : configData["MC_STEP_OFFSET"].asDouble();
	MC_NSTEPS = !configData["MC_NSTEPS"] ? MC_NSTEPS : configData["MC_NSTEPS"].asInt();
	MC_NTHERMSTEPS = !configData["MC_NTHERMSTEPS"] ? MC_NTHERMSTEPS : configData["MC_NTHERMSTEPS"].asInt();
	MC_NINITIALIZATIONSTEPS = !configData["MC_NINITIALIZATIONSTEPS"] ? MC_NINITIALIZATIONSTEPS : configData["MC_NINITIALIZATIONSTEPS"].asInt();
	MC_NADDITIONALSTEPS = !configData["MC_NADDITIONALSTEPS"] ? MC_NADDITIONALSTEPS : configData["MC_NADDITIONALSTEPS"].asInt();
	MC_NADDITIONALTHERMSTEPS = !configData["MC_NADDITIONALTHERMSTEPS"] ? MC_NADDITIONALTHERMSTEPS : configData["MC_NADDITIONALTHERMSTEPS"].asInt();
	MC_NADDITIONALINITIALIZATIONSTEPS = !configData["MC_NADDITIONALINITIALIZATIONSTEPS"] ? MC_NADDITIONALINITIALIZATIONSTEPS : configData["MC_NADDITIONALINITIALIZATIONSTEPS"].asInt();
	TIMESTEP = !configData["TIMESTEP"] ? TIMESTEP : configData["TIMESTEP"].asDouble();
	TOTALTIME = !configData["TOTALTIME"] ? TOTALTIME : configData["TOTALTIME"].asDouble();
	IMAGINARY_TIME = !configData["IMAGINARY_TIME"] ? IMAGINARY_TIME : configData["IMAGINARY_TIME"].asInt();
	ODE_SOLVER_TYPE = !configData["ODE_SOLVER_TYPE"] ? ODE_SOLVER_TYPE : configData["ODE_SOLVER_TYPE"].asInt();

	Json::Value paramsR = configData["PARAMS_REAL"];
	for (unsigned int i = 0; i < paramsR.size(); i++)
	{
		PARAMS_REAL.push_back(paramsR[i].asDouble());
	}
	Json::Value paramsI = configData["PARAMS_IMAGINARY"];
	for (unsigned int i = 0; i < paramsI.size(); i++)
	{
		PARAMS_IMAGINARY.push_back(paramsI[i].asDouble());
	}
	PARAM_PHIR = !configData["PARAM_PHIR"] ? PARAM_PHIR : configData["PARAM_PHIR"].asDouble();
	PARAM_PHII = !configData["PARAM_PHII"] ? PARAM_PHII : configData["PARAM_PHII"].asDouble();
}

void WriteConfig(string fileName, vector<double>& uR, vector<double>& uI, double phiR, double phiI)
{
	ofstream configFile;
	configFile.open(OUTPUT_DIRECTORY + fileName, ios::out);
	configFile.precision(8);

	configFile << "{" << endl;
	configFile << "\t" << "\"" << "SYSTEM_TYPE" << "\"" << " : " << "\"" << SYSTEM_TYPE << "\"," << endl;
	configFile << "\t" << "\"" << "OUTPUT_DIRECTORY" << "\"" << " : " << "\"" << originalOutputDirectory << "\"," << endl;
	configFile << "\t" << "\"" << "N" << "\"" << " : " << N << "," << endl;
	configFile << "\t" << "\"" << "DIM" << "\"" << " : " << DIM << "," << endl;
	configFile << "\t" << "\"" << "N_PARAM" << "\"" << " : " << N_PARAM << "," << endl;
	configFile << "\t" << "\"" << "RHO" << "\"" << " : " << RHO << "," << endl;
	configFile << "\t" << "\"" << "RC" << "\"" << " : " << RC << "," << endl;
	configFile << "\t" << "\"" << "MC_STEP" << "\"" << " : " << MC_STEP << "," << endl;
	configFile << "\t" << "\"" << "MC_STEP_OFFSET" << "\"" << " : " << MC_STEP_OFFSET << "," << endl;
	configFile << "\t" << "\"" << "MC_NSTEPS" << "\"" << " : " << MC_NSTEPS << "," << endl;
	configFile << "\t" << "\"" << "MC_NTHERMSTEPS" << "\"" << " : " << MC_NTHERMSTEPS << "," << endl;
	configFile << "\t" << "\"" << "MC_NINITIALIZATIONSTEPS" << "\"" << " : " << MC_NINITIALIZATIONSTEPS << "," << endl;
	configFile << "\t" << "\"" << "MC_NADDITIONALSTEPS" << "\"" << " : " << MC_NADDITIONALSTEPS << "," << endl;
	configFile << "\t" << "\"" << "MC_NADDITIONALTHERMSTEPS" << "\"" << " : " << MC_NADDITIONALTHERMSTEPS << "," << endl;
	configFile << "\t" << "\"" << "MC_NADDITIONALINITIALIZATIONSTEPS" << "\"" << " : " << MC_NADDITIONALINITIALIZATIONSTEPS << "," << endl;
	configFile << "\t" << "\"" << "TIMESTEP" << "\"" << " : " << TIMESTEP << "," << endl;
	configFile << "\t" << "\"" << "TOTALTIME" << "\"" << " : " << TOTALTIME << "," << endl;
	configFile << "\t" << "\"" << "IMAGINARY_TIME" << "\"" << " : " << IMAGINARY_TIME << "," << endl;
	configFile << "\t" << "\"" << "ODE_SOLVER_TYPE" << "\"" << " : " << ODE_SOLVER_TYPE << "," << endl;

	configFile << "\t" << "\"" << "PARAMS_REAL" << "\"" << " : " << endl;
	configFile << "\t" << "[" << endl;
	configFile << "\t\t";
	for (int i = 0; i < N_PARAM; i++)
	{
		configFile << fixed << uR[i];
		if (i != N_PARAM - 1)
		{
			configFile << ", ";
		}
	}
	configFile << endl;
	configFile << "\t" << "]," << endl;

	configFile << "\t" << "\"" << "PARAMS_IMAGINARY" << "\"" << " : " << endl;
	configFile << "\t" << "[" << endl;
	configFile << "\t\t";
	for (int i = 0; i < N_PARAM; i++)
	{
		configFile << fixed << uI[i];
		if (i != N_PARAM - 1)
		{
			configFile << ", ";
		}
	}
	configFile << endl;
	configFile << "\t" << "]," << endl;

	configFile << "\t" << "\"" << "PARAM_PHIR" << "\"" << " : " << phiR << "," << endl;
	configFile << "\t" << "\"" << "PARAM_PHII" << "\"" << " : " << phiI << endl;

	configFile << "}";

	configFile.close();
}

void WriteConfigJson(string fileName, double* uR, double* uI, double phiR, double phiI)
{
	ofstream configFile;
	Json::Value event;
	Json::Value paramsR(Json::arrayValue);
	Json::Value paramsI(Json::arrayValue);

	event["SYSTEM_TYPE"] = SYSTEM_TYPE;
	event["OUTPUT_DIRECTORY"] = OUTPUT_DIRECTORY;
	event["N"] = N;
	event["DIM"] = DIM;
	event["N_PARAM"] = N_PARAM;
	event["RHO"] = RHO;
	event["RC"] = RC;
	event["MC_STEP"] = MC_STEP;
	event["MC_STEP_OFFSET"] = MC_STEP_OFFSET;
	event["MC_NSTEPS"] = MC_NSTEPS;
	event["MC_NTHERMSTEPS"] = MC_NTHERMSTEPS;
	event["MC_NINITIALIZATIONSTEPS"] = MC_NINITIALIZATIONSTEPS;
	event["MC_NADDITIONALSTEPS"] = MC_NADDITIONALSTEPS;
	event["MC_NADDITIONALTHERMSTEPS"] = MC_NADDITIONALTHERMSTEPS;
	event["MC_NADDITIONALINITIALIZATIONSTEPS"] = MC_NADDITIONALINITIALIZATIONSTEPS;
	event["TIMESTEP"] = TIMESTEP;
	event["TOTALTIME"] = TOTALTIME;
	event["IMAGINARY_TIME"] = IMAGINARY_TIME;
	event["ODE_SOLVER_TYPE"] = ODE_SOLVER_TYPE;

	for (int i = 0; i < N_PARAM; i++)
	{
		paramsR.append(Json::Value(uR[i]));
		paramsI.append(Json::Value(uI[i]));
	}
	event["PARAMS_REAL"] = paramsR;
	event["PARAMS_IMAGINARY"] = paramsI;

	event["PARAM_PHIR"] = phiR;
	event["PARAM_PHII"] = phiI;

	configFile.open(OUTPUT_DIRECTORY + fileName, ios::out);
	configFile << event << endl;
	configFile.close();
}