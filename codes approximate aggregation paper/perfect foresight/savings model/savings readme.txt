----------------------------------------------------------------------------------------------
Approximate Aggregation Theory for Heterogeneous-agent Models
by Guillermo Hausmann-Guil (2022))
Savings model README
----------------------------------------------------------------------------------------------


This folder contains the Matlab files that replicate the solution to the Savings model of Section 6.
The files are (by type):
- Scripts:
savings_solution.m %master file that computes the perfect-foresight solution


- Functions to approximate AR1 process (not mine)
rouwen.m %Grid and transition matrix for idiosyncratic innovation

-Functions to implement the endogenous grid method
c_poli_savings_update.m %updates the consumption rule for the perfect foresight case

-Functions to compute the wealth distribution (code inside borrowed from Guerrieri and Lorenzoni (2017))
ss_distribution.m %computes the stationary distribution
update_distribution.m %updates a given current distribution  

-Functions containing systems of equations
my_savings_ss.m %to solve the steady-state equilibrium of the perfect foresight case
savings_transition.m %to compute the perfect-foresight transition for the price
variance_transition.m %to compute the perfect-foresight transition for the variance
savings_M1_Z.m %solves the M1 approximation
savings_M2_S_fun.m %to solve for the coefficients of the variance in the M2 approximation
savings_M3_S_fun.m %to solve for the coefficients of the variance and third moment in the M3 approximation
savings_M2_Z.m %solves the coefficients of the borrowing limit in the M2 approximation
savings_M3_Z.m %solves the coefficients of the borrowing limit in the M3 approximation

-Functions computing derivatives
savings_idgdX.m %computes average partial derivatives w.r.t. the sequence of borrowing limits and prices, for M1, M2, and M3 solutions 
savings_idgds.m %computes average partial derivatives w.r.t. today's savings, for M2 and M2 solutions
