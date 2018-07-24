#pragma once

#include "Utils.h"

#include <string>
#include <vector>

#include "../resources/json/json-forwards.h"
#include "../resources/json/json.h"

using namespace std;

enum ConfigItemType
{
	STRING, INT, DOUBLE, ARR_DOUBLE
};

class ConfigItem
{
public:
	string name;

	string* variableString;
	int* variableInt;
	double* variableDouble;
	vector<double>* variableArrDouble;

	ConfigItemType type;

public:
	ConfigItem(string name, string* variable, ConfigItemType type);
	ConfigItem(string name, int* variable, ConfigItemType type);
	ConfigItem(string name, double* variable, ConfigItemType type);
	ConfigItem(string name, vector<double>& variable, ConfigItemType type);

	void setValue(Json::Value value);

	string getJsonString();
};
