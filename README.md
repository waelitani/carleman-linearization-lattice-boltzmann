# Carleman Linearization of lattice Boltzmann
This is a repository of MATLAB codes accompanying "Analysis of Carleman Linearization of Lattice Boltzmann". If you find this work useful in your research, please cite:

Itani, W. & Succi, S. (2021). *Analysis of Carleman Linearization of Lattice Boltzmann*.

## LB Code
This folder contains the variants of the code used for studying and visualizing the Carleman linearization of the lattice Boltzmann scheme.

### Collision
1. carlemannBoltzmann.m: It sets up the linearized lattice Boltzmann problem with number of discrete velocites *Q*, number of dimensions *D*, lattice weights ***w***, lattice speed *c*, incompressible density *rho*, lattice vectors ***e***,  relaxation time *T*, number of timesteps *time*, ensemble size for the study *K_max*, Carleman linearization order *order*, initial flow velocity *uu*, initial discrete density values ***f_0***, timestep size *dt*. It then calls the latticeBoltzmann, carlemannLinearize, timeStep, computeResultsArray and plotSolReg functions in this order.
2. latticeBoltzmann.m: It returns the collision term ***omega*** for a single phase fluid with BGK equilibrium function for the parameters listed in carlemannBoltzmann.
3. carlemannLinearize.m: It returns the Carlemann variables *V* and matrices *A* for the different linearization orders of a given function. They are calculated recursively in increasing order to avoid computing terms that do not appear in the equations.
5. timeStep.m: It returns the function and linearized matrix to compute the next value in time for a first order equation of ***f*** driven by ***omega***. Both implicit and explicit Euler schemes are available.
6. computeResultsArray.m: It returns cell arrays containing the values of the original variables ***f*** and Carlemann variables ***V** at different times for given initial conditions, relaxation times, and linearization orders.
7. plotSolReg.m: It plots the exact and Carleman-linearized solutions of each discrete density *f_i* with the respective normalized error.

### postProcessing
The folder contains different post-processing tools to visualize the sparsity of the Carleman matrix, and calculate the norms used in the nonlinearity ratio on which the analytically-derived Carleman linearization error depends. The nonlinear driving functions ***omega1*** and ***omega2*** are assumed to pre-calculated.

### Random Initial
The folder contains a variation of the code in ***Collision*** with a post-processing function tool in dataTable.m to visualize the root sum of squares of the error with respect to the ratio of timestep to relaxation time, linearization order and distance of the initial condition to the weights vector.

## Logistic Code
The folder contains a variant of ***LB Code\Collision*** modified for a single population, namely the logistic equation with a single variable. Unlike the lattice Boltzmann carlemannLinearize.m, the function in this folder precalculates the Carlemann variables by taking the different orders of the original variable, and then computes the Carlemann linearization matrices for all chosen orders.
