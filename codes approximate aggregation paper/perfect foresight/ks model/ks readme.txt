----------------------------------------------------------------------------------------------
Approximate Aggregation Theory for Heterogeneous-agent Models
by Guillermo Hausmann-Guil (2022))
Krusell-Smith model README
----------------------------------------------------------------------------------------------


This folder contains the Matlab files that replicate the solution to the Krusell-Smith model of Section 6.
The files are (by type):
- Scripts:
ks_solution.m %master file that computes the perfect-foresight solution

-Functions to implement the endogenous grid method
c_poli_ks_update.m %updates the consumption rule

-Functions to compute the wealth distribution (code inside borrowed from Guerrieri and Lorenzoni (2017))
ss_distribution.m %computes the stationary distribution
update_distribution.m %updates a given current distribution  


-Functions containing systems of equations
my_ks_ss.m %to solve the steady-state equilibrium
my_ks_transition.m %to compute the perfect-foresight transition
ks_M1_S_fun.m %to solve the l_K coefficient from the law of motion for capital (root finding) 
ks_M1_S_fun_it.m %to solve the l_K coefficient from the law of motion for capital (iteration) 
ks_M1_Z.m %solves the l_Z coefficient from the law of motion for capital 

-Functions computing derivatives
ks_idgdX.m %computes average partial derivatives w.r.t. the sequence of capital and prices