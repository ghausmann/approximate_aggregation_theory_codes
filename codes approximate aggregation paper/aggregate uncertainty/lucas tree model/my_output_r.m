%This function updates the asset rule and the wealth distribution for a
%given current price, and a given level of output, using the
%approximate solution to form expectations about future prices.
%Inputs are P (vector of parameters), a_grid and zi (grids for asset and
%productivity values with lengths n_a and n_z), Pr (Probability transition
%matrix with size n_z-by-n_z), Yg (output grid of size n_Y), PrY
%(Probability transition matrix of size n_Y-by-n_Y), c_pol_r (consumption
%rule, a n_z-by-n_a-by-n_Y array), L (the pair (p_ss,lY)), Dt (the current
%asset distribution), ind_Yg (index for current level of Yg), x (the
%current asset price), and order (the order of approximation).

function [a_poli,Dt1] = my_output_r(P,a_grid,zi,Pr,Yg,PrY,c_pol_r,L,Dt,ind_Yg,x,order)

%Model parameters
alpha = P(3);
n_a = length(a_grid);

pt = x; %current asset price

%Compute the current consumption rule (check the function
%"c_poli_system_lucas_z_check.m" for details).
c_poli =  c_poli_lucas_r_check(P,L,a_grid,zi,Pr,Yg,ind_Yg,PrY,c_pol_r,pt,order);

%Return the corresponding rule for savings
a_poli = ((1+zi')*(1-alpha)*(Yg(ind_Yg)) + (1+a_grid)*(pt+alpha*(Yg(ind_Yg))) - c_poli)/pt - 1;
a_poli(a_poli<=-1)=-1;
a_poli(a_poli>=a_grid(n_a))=a_grid(n_a);

%Return the updated wealth distribution, using the function
%"update_distribution.m".
Dt1 = update_distribution(a_grid,Pr,a_poli,Dt);