----------------------------------------------------------------------------------------------
Approximate Aggregation Theory for Heterogeneous-agent Models
by Guillermo Hausmann-Guil (2022))
Two-asset model README
----------------------------------------------------------------------------------------------


This folder contains the Matlab codes that replicate the solution to the two-asset model of Section 7.
The files are:

local_2assets_derivatives.m %script that derives the 10 equations coming from euler-equation conditions, using symbolic differentiation
local_2assets_system.m %function that builds the 12x12 system of equations to solve the local approximation
local_2assets_solution.m %master file that computes all the results
local_2assets_errors.m %function that computes euler-equation errors
myxt.mat %Data file to replicate the simulation


