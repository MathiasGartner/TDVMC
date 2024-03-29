/*-----------------------------------------------------------------------------
 *
 * 		Name:			ConfigItem.h
 * 		Author:			Mathias Gartner
 * 		Description:	Holds a simulation parameter of type <ConfigItemType>.
 * 						Simulations parameters are typically read from a
 * 						configuration file in JSON format.
 *
 * --------------------------------------------------------------------------*/

#pragma once

#include "Utils.h"

#include <string>
#include <vector>

#include "../resources/json/json-forwards.h"
#include "../resources/json/json.h"

using namespace std;

enum ConfigItemType
{
	STRING, INT, DOUBLE, ARR_INT, ARR_DOUBLE
};

class ConfigItem
{
public:
	string name;
	bool allowChangeAtRuntime;

	string* variableString;
	int* variableInt;
	double* variableDouble;
	vector<int>* variableArrInt;
	vector<double>* variableArrDouble;

	ConfigItemType type;

private:
	ConfigItem(string name, ConfigItemType type, bool allowChangeAtRuntime = false);

public:
	ConfigItem(string name, string* variable, ConfigItemType type, bool allowChangeAtRuntime = false);
	ConfigItem(string name, int* variable, ConfigItemType type, bool allowChangeAtRuntime = false);
	ConfigItem(string name, double* variable, ConfigItemType type, bool allowChangeAtRuntime = false);
	ConfigItem(string name, vector<int>& variable, ConfigItemType type, bool allowChangeAtRuntime = false);
	ConfigItem(string name, vector<double>& variable, ConfigItemType type, bool allowChangeAtRuntime = false);

	void setValue(Json::Value value);

	string getJsonString();
};
