%This function uses the endogenous grid method, mapping an old consumption
%rule to a new one. When both are the same, and the asset price is constant
%over time, we have the optimal SS rule. The function can also be used to
%perform "backward iteration", which allows to calculate approximate
%derivatives of today's policy rules with respect to future prices and
%shocks.
%Inputs are P (vector of parameters), a_grid and zi (grids for asset of
%length n_a and productivity values of length n_z), Pr (Probability transition
%matrix of size n_z-by-n_z), c_pol (the old/tomorrow's consumption rule, a n_z-by-n_a
%matrix), Y0 and Y1 (today's and tomorrow's output), and p0 and
%p1 (today's and tomorrow's asset prices).
%When computing the SS optimal rule, simply set Y0=Y1=1, and p0=p1=p_ss.

function y = c_poli_lucas_update(P,a_grid,zi,Pr,c_pol,Y0,Y1,p0,p1)

%Model parameters
beta = P(1);
gamma = P(2);
alpha = P(3);

%lenghts of the asset and productivity grids
n_a = length(a_grid);
n_z = length(zi);

c_poli = zeros(n_z,n_a);

%Compute today's consumption using the Euler equation, for fixed values of
%savings, prices, and future excess output.
c_i = ((1/(  (beta/(p0)  )*( (p1) + alpha*(Y1) )))./(  Pr*((1./c_pol).^gamma)   )  ).^(1/gamma);


%Compute today's implied initial level of assets.
a_grid1 = ( c_i + (p0)*(1+a_grid) - (1+zi')*(1-alpha)*(Y0) )/((p0) + alpha*(Y0))  -1;


%Interpolate to complete the mapping from the asset grid to the new
%consumption rule
for j=1:n_z
    
    c_poli(j,:) = interp1(a_grid1(j,:), c_i(j,:), a_grid, 'linear', 'extrap');
        
end

%Correct values that violate the constraints (-1<=a_poli<=a_grid(I)):

%Compute the savings rule
a_poli = ( (1+zi')*(1-alpha)*(Y0) + (1+a_grid)*((p0)+alpha*(Y0)) - c_poli)/p0 - 1;
%Consumption if the borrowing constraint binds
c_poli_cd = (1+zi')*(1-alpha)*(Y0) + ((p0)+alpha*(Y0))*(1+a_grid);
%Consumption if the upper limit for savings binds
c_poli_cu = (1+zi')*(1-alpha)*(Y0) + ((p0)+alpha*(Y0))*(1+a_grid) - (p0)*(1+a_grid(n_a));

%Replace values that violate constraints by the proper ones
c_poli(a_poli<=-1) = c_poli_cd(a_poli<=-1);
c_poli(a_poli>=a_grid(n_a)) = c_poli_cu(a_poli>=a_grid(n_a));

%Return today's policy rule
y = c_poli;