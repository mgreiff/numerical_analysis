# numerical_analysis
This repository contains a collection of useful algorithms in everything from
optimization to parallell PDE solvers. Each subdirectory has one or many
complete reports located in the /report directories, or a simple README.md file
if the extent of the project doesn't warrant a full report.

## Contents
#### /optimization_methods
* A solver implementing quasi-newton methods (BFGS/DFP) for unconstrained
  optimization (similar to matlabs fminunc but without the need for costly licenses).
  The repo is complete with three examples of barrier/penalty formulations of
  constrained convex problems and a detailed report. (**Matlab**)
* Solvers for the travelling salesman problem (TSP), both with a genetic algorithm
  (not optimal but fast) and with a branch and bound method (optimal) with nice
  visualization options. (**Python**)
* A solver for the resource allocation problem (RAP), using a dual shadow cost
  formulation. (**Python**)

#### /multigrid_methods
* A multigrid helmholtz equation solver complete with an example simulating the
  two dimensional wave equation with a vibrating membrane on the unit square, complete
  with a detailed report. (**Matlab**)

#### /parallell_programming_methods
* A simple parallell solver for heat propagation in a irregular geometry, here a
  two dimensional ''house'' with three rooms. The scipt shows how the Pythons MPI
  extension can be used to solve the Poisson equation in a tricky computational
  domain. (**Python**)
