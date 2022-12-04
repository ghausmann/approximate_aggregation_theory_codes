%This function computes the average derivative of the asset rule w.r.t
%output when the later is given by a grid, using the central difference
%formula.
%Inputs are P (vector of parameters), a_grid, zi, and DYg (grids for asset
%of length n_a, productivity values of length n_z, and (Y-1) values of length
%n_Y), c_pol_r (the consumption rule, a n_z-by-n_a-by-n_Y array), Dss (the
%stationary distribution), dS (the differentiation step), and x (the vector
%of coefficients pss and lY).

function y = lucas_M1_r_dgdY_fo(P,a_grid,zi,DYg,c_pol_r,Dss,dS,x)

%Model parameters
alpha = P(3);

%Asset price coefficients
pss = x(1);
lY = x(2);

dSd = -dS; %the small negative increase

%lenghts of the asset, productivity, and output grids
n_a = length(a_grid);
n_z = length(zi);
n_Y = length(DYg);

%Interpolate along the z dimension, to approximate the consumption rule
%when z changes by a small amount (positive and negative)
c_pert = permute(c_pol_r,[3 2 1]);
c_pert_rep = reshape(c_pert,n_Y, n_z*n_a);

%Compute the asset rule for a small increase in Y
c_pert_rep_dSu = interp1(DYg,c_pert_rep,dS,'spline'); %positive increase
c_poli_dSu = (reshape(c_pert_rep_dSu',n_a,n_z))'; %new consumption rule
a_poli_dSu = ((1+zi')*(1-alpha)*(1+dS) + (1+a_grid)*((pss + lY*dS )+alpha*(1+dS)) - c_poli_dSu)/(pss + lY*dS ) - 1;

%Compute the asset rule for a small decrease in Y
c_pert_rep_dDd = interp1(DYg,c_pert_rep,dSd,'spline'); %negative increase
c_poli_dSd = (reshape(c_pert_rep_dDd',n_a,n_z))'; %new consumption rule
a_poli_dSd = ((1+zi')*(1-alpha)*(1+dSd) + (1+a_grid)*((pss + lY*dSd )+alpha*(1+dSd)) - c_poli_dSd)/(pss + lY*dSd ) - 1;

%Approximate the derivative using the central difference formula
dgdY = (0.5*(a_poli_dSu-a_poli_dSd)/dS);

%Return the average derivative
y = sum(sum(Dss.*dgdY));