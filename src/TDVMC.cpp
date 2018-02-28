/*
 ============================================================================
 Name        : TDVMC.c
 Author      : Mathias Gartner
 Version     :
 Copyright   : 
 Description :
 ============================================================================
 */

#include "mpi.h"
#include <iostream>
#include "test/PiCalculator.h"
#include <stdlib.h>

using namespace std;

int main(int argc, char **argv) {
	int val = -1;
	cout << "start" << endl;

	if (argc > 1)
	{
		if (strcmp(argv[1], "-pitest") == 0)
		{
			CalculatePi(argc, argv);
		}
		else
		{
			//TODO: TDVMC
		}
	}

	val = 0;
	cout << "finish " << val << endl;
	return val;
}

