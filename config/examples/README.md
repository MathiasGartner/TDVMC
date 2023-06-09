# Description of fields in config files

- the config file is given in JSON notation
- for boolean types the values 0 (=false) and 1 (=true) are used.
- lists are specified with [ ... ]

## the following properties need to be set in the config file:

### general properties

- `CONFIG_VERSION` : the currently supported version of a config file is specified in src/TDVMC.cpp in the variable `requiredConfigVersion` and is currently set to "0.24". If new config fields are added/renamed/deleted, this value should be updated.
- `SYSTEM_TYPE` : the type of physical model that is simulated. The value must correspond to a class located in src/PhysicalSystems/ eg. "Bosons1D"
- `OUTPUT_DIRECTORY` : an **existing** directory in which a new simulation directory for all the output files will be created. The simulation directory is created as "step=<MC_NSTEPS>_therm=<MC_NTHERMSTEPS>_time=<TIMESTEP><OUT_DIR_SUFFIX>" if no value for <OUT_DIR_NAME> is specified.
- `OUT_DIR_SUFFIX` : name that is appended to the created simulation directory, eg: "_GS_1D_N=100_"
- `OUT_DIR_NAME` : if specified the created simulation directory will be named exactly by this eg. "1D_groundstate_N=40"

### realted to simulated system, wavefunction

- `N` : number of particles, eg. 100
- `LBOX` : size of the simulation box, eg. 100.0
- `DIM` : number of dimensions, eg. one of {1, 2, 3}
- `N_PARAM` : number of variational parameters, eg. 100
- `USE_PARAM_START` : the index at which the parameter array will be used. mainly used for code development (eg. 10 if the first 10 parameters should stay at a fixed value). default value. 0
- `USE_PARAM_END` : same functionality as USE_PARAM_START. mainly used for code development. default value is N_PARAM - 1. eg 99
- `RHO` : the particle density. This value is ignored in most SYSTEM_TYPE implementations since the density is given by N/LBOX
- `RC` : not used anymore. has been used to specify the position r_c where the Lennard-Jones potential is cut off

### Monte Carlo parameters

- `MC_STEP` : magnitude of coordinate displacement during a Monte Carlo step
- `MC_STEP_OFFSET` : offset of coordinate displacement for a Monte Carlo step
- `MC_NSTEPS` : Monte Carlo samples for a timestep created in each process
- `MC_NTHERMSTEPS` : intermediate MC samples generated between each of the MC_NSTEPS steps
- `MC_NINITIALIZATIONSTEPS` : initial MC samples generated before the first of the MC_NSTEPS
- `MC_VERY_FIRST_NINITIALIZATIONSTEPS` : initial MC samples generated before the very first timestep of the simulation
- `MC_NADDITIONALSTEPS` : every CALCULATE_ADDITIONAL_DATA_EVERY_NTH_STEP timesteps additional data (eg. density, pair distribution, static structure factor, ...) are calculated. for this MC_NADDITIONALSTEPS Monte Carlo samples are used
- `MC_NADDITIONALTHERMSTEPS` : thermalization steps for MC_NADDITIONALSTEPS
- `MC_NADDITIONALINITIALIZATIONSTEPS` : initialization steps for MC_NADDITIONALSTEPS
- `MC_NFINALSTEPS_MULTIPLICATOR` : at the very end of the simulation the additional observables are estimated a last time. the Monte Carlo samples used for this last evaluation can be increased by the factor MC_NFINALSTEPS_MULTIPLICATOR in higher accuracy is needed

### time propagation

- `TIMESTEP` : the timestep for time propagation
- `TOTALTIME` : the total simulation time
- `IMAGINARY_TIME` : whether the system is propagated in imaginary or real time (0 = imaginary, 1 = real)
- `ODE_SOLVER_TYPE` : type of ODE solver for time propagation. commonly used values are 0 (simple Euler method for imaginary time) or 2 (RungeKutta4 for real time evolution)
- `LINEAR_EQUATION_SOLVER_TYPE` : 0 (Cholesky), 1 (SVD) 
- `USE_PRECONDITIONING` : use preconditioning of the equation system. 0 (false), 1 (true)
- `USE_PARAMETER_ACCEPTANCE_CHECK` : check if the parameter update is meaningful. 0 (false), 1 (true)
- `PARAMETER_ACCEPTANCE_CHECK_TYPE` : which type of parameter acceptance check should be used. 0 (None), 1 (if parameters are finite), 2 (relative change beyond threshold), 3 (relative change of the last 10 steps within threshold), 4 (energy and parameter below threshold), 5 (parameters finite, energy difference below threshold)

### observables during time propagation

- `WRITE_EVERY_NTH_STEP_TO_FILE` : write data every WRITE_EVERY_NTH_STEP_TO_FILE to files
- `CALCULATE_ADDITIONAL_DATA_EVERY_NTH_STEP` : calculate additional proeprties like density, pair distribution, static structure factor, ... every CALCULATE_ADDITIONAL_DATA_EVERY_NTH_STEP steps. should be the same as WRITE_EVERY_NTH_STEP_TO_FILE
- `WRITE_SINGLE_FILES` : write a single file for every timestep with current timestep append to filenames. 0 (no), 1 (yes)
- `MC_NSTEP_MULTIPLICATION_FACTOR_FOR_WRITE_DATA` : if data is written to file, increase the used Monte Carlo samples by this factor to get better accuracy on this written data. default 1

### numerics related

- `USE_MEAN_FOR_FINAL_PARAMETERS` : used for imaginary time simulations (ground state). the final parameters that are saved to file are the mean over the last timesteps. be careful with this - this option can lead to not very good ground state wavefunctions. default: 0
- `USE_NORMALIZE_WF` : if the wavefunction is normalized in every time step
- `USE_ADJUST_PARAMETERS` : if the parameters should be adjusted in every step (eg. to have mean 0). default 0 (no adjustment)
- `UPDATE_SAMPLES_EVERY_NTH_STEP` : specify a value larger than 0 if generated Monte Carlo samples should be used for successive timesteps. If a value of 0 is specified no samples are reused. This is the standard behaviour. Values other than 0 are useful if generating single samples is very expensive.
- `UPDATE_SAMPLES_PERCENT` : how much (in %) of the reused samples should be updated in every timestep

### output data related

- `GR_BIN_COUNT` : bins for creating a histogram in the evaulation of the pair distributin function g_2(r_ij)
- `RHO_BIN_COUNT` : bins for creating a histogram in the evaulation of the density rho(r)

### non uniform grid for splines

- `USE_NURBS` : use Non-Uniform-Rational-B-Splines: 0 (no), 1 (yes). default is no
- `NURBS_GRID` : a list of the NURBS grid points

### particle mixtures

- `PARTICLE_TYPES` : for systems with different particle types a list with up to N entries (integer values) can be specified. If less than N values are given, the last value is repeated for the remaining particles. The specific implementation for each particle type is according to the SYSTEM_TYPE. eg [1, 2, 0] for one particle of type 1, one particle of type 2 and the remaining particles of type 0

### additional values specific to the simulated system

- `SYSTEM_PARAMS` : a list of SYSTEM_TYPE specific parameters. Mostly used for interaction parameters or external potential parameters. eg. [1.0, 10.0, 1.0, 2.0]

### initial values for variational parameters

- `PARAMS_REAL` : a list for the initial values of the real part of the variational parameters. If [ 0.0 ] is specified, all parameters are initialized to zero.
- `PARAMS_IMAGINARY` : a list for the initial values of the imaginary part of the variational parameters. If [ 0.0 ] is specified, all parameters are initialized to zero.
- `PARAM_PHIR` : initial value for the normalization parameter Re[\phi]
- `PARAM_PHII` : initial value for the phase parameter Im[\phi]
