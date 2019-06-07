//!HFD-B(He) Aziz et al: Mol. Phys. 61, 1487 (1987)
#pragma once

#include "IPairPotential.h"

namespace Potentials
{

class KTTY_He_Na: public IPairPotential
{
public:
	double epsil;
	double d;
	double b1;
	double b2;
	double c6;
	double c8;
	double c10;

	double c12;
	double c14;
	double c16;

	double fak2;
	double fak3;
	double fak4;
	double fak5;
	double fak6;
	double fak7;
	double fak8;
	double fak9;
	double fak10;
	double fak11;
	double fak12;
	double fak13;
	double fak14;
	double fak15;
	double fak16;

public:
	KTTY_He_Na();

	double GetPotential(double r);
};

}
