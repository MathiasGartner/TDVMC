# tVMC - Time Dependent Variational Monte Carlo

source code using in: https://epub.jku.at/obvulihs/id/8568071

## building from source using Makefile
you can specify options according to the machine your are using for compiling the project, eg.
- CXX = {mpic++, mpicxx, ...}
- C++ standard >= -std=c++17
- build directory (default is ./build)

## config file for simulation run
- config files are in JSON format
- example config file can be found in directory config/examples/
- there is also a README file config/examples/README.md with explanation for all the fields in the config
- at the end of every simulation a file AAfinish_tVMC.config is generated in the simulation directory. This file contains the final parameters and can be used in a subsequent simulation (eg. for further imaginary time propagation if parameters are not converged, for real time propagation following an imaginary time propagation, for further real time simulations, ...).

## the most advanced PhysicalSystem implemetations are
- `InhContactBosons`: Bosons in 1D with contact interaction and external potential. Periodic boundary conditions. Wavefunction has single particle and pair correlation parts
- `BosonMixtureCluster`: Bosonic mixtures in 3D. Used for ground state optimization of few-body clusters. Piecewise defined wavefunction.
- `Bosons1D`: Bosons in 1D with arbitrary interaction potential. Homogeneous system with periodic boundary conditions. Wavefunction has only pair correlations.
- `Bosons1DMixture`: Bosonic mixtures in 1D with arbitrary interaction potential. Homogeneous system with periodic boundary conditions. Wavefunction has only pair correlations.

## external resources
JSON parser: https://github.com/open-source-parsers/jsoncpp
Eigen library: https://eigen.tuxfamily.org
