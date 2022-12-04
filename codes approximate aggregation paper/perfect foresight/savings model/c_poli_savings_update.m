%This function uses the endogenous grid method, mapping an old consumption
%rule to a new one. When both are the same, and the bond price is constant
%over time, we have the optimal SS rule. The function can also be used to
%perform "backward iteration".
%Inputs are P (vector of parameters), b_grid and ui (grids for bonds of
%length n_b and income shock values of length n_u), Pr (Probability
%transition matrix of size n_u-by-n_u), c_pol (the old/tomorrow's
%consumption rule, a n_u-by-n_b matrix), phi_t (current borrowing limit),
%and qt (current bond price).
% When computing the SS optimal rule, simply set phi_t=phi_ss, and qt=q_ss.

function y = c_poli_savings_update(P,b_grid,ui,Pr,c_pol,phi_t,qt)

%Parameter values
beta = P(1);
gamma = P(2);

%lenghts of the asset and income grids
n_b = length(b_grid);
n_u = length(ui);

c_poli = zeros(n_u,n_b);

%Compute today's consumption using the Euler equation, for fixed values of
%today's savings and income.
c_i = ((1/(beta/(qt)))./(     Pr*((1./c_pol).^gamma)   )     ).^(1/gamma);

%Compute today's implied initial level of assets.
b_grid1 = c_i + (qt)*b_grid - exp(ui');

%Interpolate to complete the mapping from the asset grid to the new
%consumption rule
for j=1:n_u
    
    c_poli(j,:) = interp1(b_grid1(j,:), c_i(j,:), b_grid, 'linear', 'extrap');
    
end

%Correct values that violate the constraints (0<=a_poli<=a_grid(I)):

%Compute the savings rule
b_poli = ( exp(ui') + b_grid - c_poli)/qt;
%Consumption if the borrowing constraint binds
c_poli_cd = exp(ui') + b_grid + (qt )*phi_t;
%Consumption if the upper limit for savings binds
c_poli_cu = exp(ui') + b_grid - (qt )*b_grid(n_b);

%Replace values that violate constraints by the proper ones
c_poli(b_poli<=-phi_t) = c_poli_cd(b_poli<=-phi_t);
c_poli(b_poli>=b_grid(n_b)) = c_poli_cu(b_poli>=b_grid(n_b));

%Return today's policy rule
y = c_poli;