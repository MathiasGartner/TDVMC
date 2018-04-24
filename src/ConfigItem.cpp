#include "ConfigItem.h"

ConfigItem::ConfigItem(string name, string* variable, ConfigItemType type)
{
	this->name = name;
	this->variableString = variable;
	this->type = type;
}

ConfigItem::ConfigItem(string name, int* variable, ConfigItemType type)
{
	this->name = name;
	this->variableInt = variable;
	this->type = type;
}

ConfigItem::ConfigItem(string name, double* variable, ConfigItemType type)
{
	this->name = name;
	this->variableDouble = variable;
	this->type = type;
}

ConfigItem::ConfigItem(string name, vector<double>& variable, ConfigItemType type)
{
	this->name = name;
	this->variableArrDouble = &variable;
	this->type = type;
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
	else if (this->type == ARR_DOUBLE)
	{
		s = "\n\t[\n\t\t " + JoinVector(*variableArrDouble) + "\n\t]";
	}
	return s;
}
