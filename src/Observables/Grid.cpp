/*-----------------------------------------------------------------------------
 *
 * 		Name:			Grid.cpp
 * 		Author:			Mathias Gartner
 * 		Description:	see Grid.h
 *
 * --------------------------------------------------------------------------*/

#include "Grid.h"

#include "../Utils.h"

namespace Observables
{
Grid::Grid()
{
	this->name = "";
	this->min = 0.0;
	this->max = 0.0;
	this->count = 0;
	this->spacing = 0.0;
}

void Grid::Init(double min, double max, double spacing)
{
	this->min = min;
	this->max = max;
	this->spacing = spacing;
	this->count = (this->max - this->min) / this->spacing;
	InitVector(this->grid, this->count, 0.0);

	double val = this->min;
	for (unsigned int i = 0; i < grid.size(); i++)
	{
		grid[i] = val;
		val += this->spacing;
	}
	this->InitGridCenterPoints();
}

void Grid::Init(vector<double> g)
{
	this->min = g[0];
	this->max = g[g.size() - 1];
	this->spacing = 0;
	this->count = g.size();
	this->grid = g;
}

void Grid::InitGridCenterPoints()
{
	InitVector(this->gridCenterPoints, this->count, 0.0);
	for (int i = 0; i <= this->count - 2; i++)
	{
		this->gridCenterPoints[i] = (this->grid[i + 1] + this->grid[i]) / 2.0;
	}
	this->gridCenterPoints[this->count - 1] = (this->max + this->grid[this->count - 1]) / 2.0;
}

int Grid::GetIndex(double value)
{
	//INFO: only for grid that starts at zero. otherwise: (gridValue - this->gridMin) has to be used.
	//INFO: recalculate if spacing is not uniform
	int index = 0;
	double interval;
	interval = value / this->spacing;
	index = floor(interval);
	return index;
}

}
