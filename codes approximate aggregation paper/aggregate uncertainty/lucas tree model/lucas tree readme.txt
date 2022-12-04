----------------------------------------------------------------------------------------------
Approximate Aggregation Theory for Heterogeneous-agent Models
by Guillermo Hausmann-Guil (2022))
Lucas-tree model README
----------------------------------------------------------------------------------------------


This folder contains the Matlab files that replicate the solution to the Lucas-tree model of Section 7.
The files are (by type):
- Scripts:
lucas_pf_solution.m %computes the perfect-foresight solution
lucas_r_solution_fo.m %computes the first-order solution with aggregate risk
lucas_r_solution_so.m %computes the second-order solution with aggregate risk

- Functions to approximate AR1 process (not mine)
rouwen.m %Grid and transition matrix for aggregate output
hitm_s.m %simulates the process

-Functions to implement the endogenous grid method
c_poli_lucas_update.m %updates the consumption rule for the perfect foresight case
c_poli_lucas_r_update.m %updates the consumption rule for the aggregate risk case (first and second order)
c_poli_lucas_r_check.m %updates the consumption rule for the accuracy test

-Functions to compute the wealth distribution (code inside borrowed from Guerrieri and Lorenzoni (2017))
ss_distribution.m %computes the stationary distribution
update_distribution.m %updates a given current distribution  
my_output_r.m %updates the asset rule and current distribution for the accuracy test

-Functions containing systems of equations
my_lucas_ss.m %to solve the steady-state equilibrium of the perfect foresight case
my_lucas_pf_transition.m %to compute the perfect-foresight transition
lucas_M1_r_fo.m %to solve the first-order approximation with aggregate risk 
lucas_M1_r_so.m %to solve the second-order approximation with aggregate risk 
my_check_r.m %to solve for the current asset price clearing the market in the accuracy test
lucas_M1_pf.m %solves the M1 approximation with perfect foresight

-Functions computing derivatives
lucas_idgdX.m %computes average partial derivatives w.r.t. the sequence of ouput and prices, for M1 solution with perfect foresight
lucas_M1_r_dGdY_fo.m %computes the average (total) derivative w.r.t. current ouput, for first-order solution with aggregate risk
lucas_M1_r_dGdY_so.m %computes the average (total) derivative w.r.t. current ouput, for second-order solution with aggregate risk

-MAT files
my_hit_baseline.mat %to simulate the economy for the accuracy test
my_hit_highrisk.mat %to simulate the economy for the accuracy test

