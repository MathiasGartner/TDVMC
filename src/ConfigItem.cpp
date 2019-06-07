#include "ConfigItem.h"

ConfigItem::ConfigItem(string name, ConfigItemType type, bool allowChangeAtRuntime)
{
	this->name = name;
	this->variableString = 0;
	this->variableInt = 0;
	this->variableDouble = 0;
	this->variableArrInt = 0;
	this->variableArrDouble = 0;
	this->type = type;
	this->allowChangeAtRuntime = allowChangeAtRuntime;
}

ConfigItem::ConfigItem(string name, string* variable, ConfigItemType type, bool allowChangeAtRuntime) :
		ConfigItem(name, type, allowChangeAtRuntime)
{
	this->variableString = variable;
}

ConfigItem::ConfigItem(string name, int* variable, ConfigItemType type, bool allowChangeAtRuntime) :
		ConfigItem(name, type, allowChangeAtRuntime)
{
	this->variableInt = variable;
}

ConfigItem::ConfigItem(string name, double* variable, ConfigItemType type, bool allowChangeAtRuntime) :
		ConfigItem(name, type, allowChangeAtRuntime)
{
	this->variableDouble = variable;
}

ConfigItem::ConfigItem(string name, vector<int>& variable, ConfigItemType type, bool allowChangeAtRuntime) :
		ConfigItem(name, type, allowChangeAtRuntime)
{
	this->variableArrInt = &variable;
}

ConfigItem::ConfigItem(string name, vector<double>& variable, ConfigItemType type, bool allowChangeAtRuntime) :
		ConfigItem(name, type, allowChangeAtRuntime)
{
	this->variableArrDouble = &variable;
}

void ConfigItem::setValue(Json::Value value)
{
	if (value != 0)
	{
		if (this->type == STRING)
		{
			*variableString = value.asString();
		}
		else if (this->type == INT)
		{
			*variableInt = value.asInt();
		}
		else if (this->type == DOUBLE)
		{
			*variableDouble = value.asDouble();
		}
		else if (this->type == ARR_INT)
		{
			variableArrInt->clear();
			for (unsigned int i = 0; i < value.size(); i++)
			{
				int v = value[i].asInt();
				variableArrInt->push_back(v);
			}
		}
		else if (this->type == ARR_DOUBLE)
		{
			variableArrDouble->clear();
			for (unsigned int i = 0; i < value.size(); i++)
			{
				double d = value[i].asDouble();
				variableArrDouble->push_back(d);
			}
		}
	}
}

string ConfigItem::getJsonString()
{
	string s = "";
	if (this->type == STRING)
	{
		s = "\"" + *variableString + "\"";
	}
	else if (this->type == INT)
	{
		s = to_string(*variableInt);
	}
	else if (this->type == DOUBLE)
	{
		s = to_string(*variableDouble);
	}
	else if (this->type == ARR_INT)
	{
		s = "\n\t[\n\t\t " + JoinVector(*variableArrInt) + "\n\t]";
	}
	else if (this->type == ARR_DOUBLE)
	{
		s = "\n\t[\n\t\t " + JoinVector(*variableArrDouble) + "\n\t]";
	}
	return s;
}
