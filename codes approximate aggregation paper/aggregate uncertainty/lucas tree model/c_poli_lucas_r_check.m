%This function uses the endogenous grid method, mapping an old consumption
%rule to a new one, accounting for both individual and aggregate risk. It
%is similar to the function "c_poli_lucas_r_update.m", except that the
%current asset price is a given number, and that it returns only the
%consumption rule conditional to a given value of output (thus a n_z-by-n_a
%matrix, and not the whole array).
%Inputs are P (vector of parameters), L (the coefficients of the price
%function), a_grid and zi (grids for asset of length n_a and productivity
%values of length n_z), Pr (Probability transition matrix of size
%n_z-by-n_z), Yg (the grid for output of size n_Y), ind_Yg (the index of
%current value of Yg), PrY (Probability transition matrix of size
%n_Y-by-n_Y), c_pol (tomorrow's consumption rule, a n_z-by-n_a-by-n_Y
%array), pt (the current asset price), and order (the order of
%approximation).

function y = c_poli_lucas_r_check(P,L,a_grid,zi,Pr,Yg,ind_Yg,PrY,c_pol,pt,order)

%Model parameters
beta = P(1);
gamma = P(2);
alpha = P(3);

%Price function coefficients
p_ss = L(1);
lY = L(2);

%Include the quadratic coefficient unless order==1
if order==1
    lY2 = 0;
else
    lY2 = L(3);
end
    

%lenghts of the asset, productivity, and output grids
n_a = length(a_grid);
n_z = length(zi);
n_Y = length(Yg);

c_poli = zeros(n_z,n_a);

%select the proper conditional distribution for tomorrow's 
%output
PrY_cond = PrY(ind_Yg,:);
PrY_cond_res = reshape(PrY_cond,1,1,n_Y);

%array for conditional expected marginal utilities
Prc_pol = zeros(n_z,n_a,n_Y);

%Compute each expected marginal utility, conditional on Y
for lly=1:n_Y
    
    Prc_pol(:,:,lly) = ((p_ss + lY*(Yg(lly)-1) + 0.5*lY2*(Yg(lly)-1)^2 ) + alpha*( Yg(lly) ))*Pr*(1./c_pol(:,:,lly)).^gamma;
    
end

%Compute the unconditional expected marginal utility
sPrc_cpol = sum(PrY_cond_res.*Prc_pol,3);

%Compute today's consumption using the Euler equation, for fixed values of
%today's savings and the current asset price pt.
c_i = ((1/(  (beta/(pt )  )))./(     sPrc_cpol   )     ).^(1/gamma);

%Compute today's implied initial level of assets.
a_grid1 = ( c_i + (pt)*(1+a_grid) - (1+zi')*(1-alpha)*(Yg(ind_Yg)) )/((pt) + alpha*(Yg(ind_Yg))) -1;

%Interpolate to complete the mapping from the asset grid to the new
%consumption rule
for j=1:n_z
    
    c_poli(j,:) = interp1(a_grid1(j,:), c_i(j,:), a_grid, 'linear', 'extrap');
    
end

%Correct values that violate the constraints (-1<=a_poli<=a_grid(I)):

%Compute the savings rule
a_poli = (  (1+zi')*(1-alpha)*(Yg(ind_Yg)) + ((pt)+alpha*(Yg(ind_Yg)))*(1+a_grid)   -  c_poli )/(pt) - 1;
%Consumption if the borrowing constraint binds
c_poli_d = (1+zi')*(1-alpha)*(Yg(ind_Yg)) + ((pt)+alpha*(Yg(ind_Yg)))*(1+a_grid) ;
%Consumption if the upper limit for savings binds
c_poli_u = (1+zi')*(1-alpha)*(Yg(ind_Yg)) + ((pt)+alpha*(Yg(ind_Yg)))*(1+a_grid) - (pt)*(1+a_grid(n_a));

%Replace values that violate constraints by the proper ones
c_poli(a_poli<=-1) = c_poli_d(a_poli<=-1);
c_poli(a_poli>=a_grid(n_a)) = c_poli_u(a_poli>=a_grid(n_a));

%Return today's policy rule, for today's value of Yg
y = c_poli;