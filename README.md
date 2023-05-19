# tVMC - Time Dependent Variational Monte Carlo

Source code used in:

- Mathias Gartner, Ferran Mazzanti, Robert E. Zillich  
  "Time-dependent variational Monte Carlo study of the dynamic response of bosons in an optical lattice"  
  SciPost Phys. 13, 025 (2022)  
  https://scipost.org/10.21468/SciPostPhys.13.2.025
- Mathias Gartner, PhD thesis  
  "Monte Carlo simulations of non-equilibrium dynamics in bosonic many-body systems"  
  https://epub.jku.at/obvulihs/id/8568071
- Mathias Gartner, David Miesbauer, Michael Kobler, Julia Freund, Giuseppe Carleo, Robert E. Zillich  
  "Interaction quenches in Bose gases studied with a time-dependent hypernetted-chain Euler-Lagrange method"  
  https://arxiv.org/abs/2212.07113

## Building from source using Makefile

This program uses MPI for multiprocess execution. You may need to specify options int the `Makefile` according to the machine your are using for compilation, eg.

- CXX = {mpic++, mpicxx, ...}
- C++ standard >= -std=c++17
- build directory (default is ./build)

## Run a simulation

The simulation config file needs to be provided as the first parameter to the program.

`./build/TDVMC ./config/Bosons1D_quench.config`

Start the program with multiple processes using the required MPI command on the host machine, eg:

`mpirun -np 250 ./build/TDVMC ./config/Bosons1D_quench.config`

`mpiexec -perhost 250 ./build/TDVMC ./config/Bosons1D_quench.config`

## Configuration file for simulation run

- config files are in JSON format
- example config file can be found in directory `config/examples/`
- there is also a README file `config/examples/README.md` with explanation for all the fields in the config
- at the end of every simulation a file `AAfinish_tVMC.config` is generated in the simulation directory. This file contains the final parameters and can be used in a subsequent simulation (eg. for further imaginary time propagation if parameters are not converged, for real time propagation following an imaginary time propagation, for further real time simulations, ...).

## The most advanced PhysicalSystem implemetations are

- `InhContactBosons`: Bosons in 1D with contact interaction and external potential. Periodic boundary conditions. Wavefunction has single particle and pair correlation parts
- `BosonMixtureCluster`: Bosonic mixtures in 3D. Used for ground state optimization of few-body clusters. Piecewise defined wavefunction.
- `Bosons1D`: Bosons in 1D with arbitrary interaction potential. Homogeneous system with periodic boundary conditions. Wavefunction has only pair correlations.
- `Bosons1DMixture`: Bosonic mixtures in 1D with arbitrary interaction potential. Homogeneous system with periodic boundary conditions. Wavefunction has only pair correlations.

## Adding new types of PhysicalSystems

To simulate new types of physical systems a new class has to be added in `src/PhysicalSystems` which implements the interface `src/PhysicalSystems/IPhysicalSystem.h`. This new class then needs to be registed in `src/TDVMC.cpp` with the following steps:

- include the class via `#include "PhysicalSystems/NewSystem.h"`
- add initialization code in method `InitializePhysicalSystem()`
- add file output specifications in method `WriteGridFiles()` and `WriteAdditionalObservablesToFiles()`
- optional: add code for local testing in `main()`

## External resources

- JSON parser: https://github.com/open-source-parsers/jsoncpp
- Eigen library: https://eigen.tuxfamily.org
