# tVMC
tVMC - Time Dependent Variational Monte Carlo

## resources
use json parser from https://github.com/open-source-parsers/jsoncpp

## building from source using Makefile
you can specify options according to the machine your are compiling the project, eg.
- CXX = {mpic++, mpicxx, ...}
- C++ standard >= -std=c++17
- build directory (default is ./build)

## config file for simulation run
- config files are in JSON format
- example config file can be found in directory config/examples/
- there is also a README file config/examples/README.md with explanation for all the fields in the config
- at the end of every simulation a file AAfinish_tVMC.config is generated in the simulation directory. This file contains the final parameters an can be used in a subsequent simulation (eg. for further imaganinary time propagation if parameters are not converged, for real time propagation following an imaginary time propagation, for further real time simulations, ...). 

