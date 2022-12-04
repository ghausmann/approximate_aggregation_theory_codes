%This function uses the endogenous grid method, mapping an old consumption
%rule to a new one. When both are the same, and the capital rental price is
%constant over time, we have the optimal SS rule. The function can also be
%used to perform "backward iteration" of the policy rule.
%Inputs are P (vector of parameters), a_grid and zi (grids for asset shares
%of length n_a and productivity values of length n_z), Pr (Probability
%transition matrix of size n_z-by-n_z), c_pol (the old/tomorrow's consumption
%rule, a n_z-by-n_a matrix), r0 and r1 (today's and tomorrow's rental
%prices), w0 (today's wage rate) and K0 and K1 (today's and tomorrow's
%aggregate capital).
%When computing the SS optimal rule, simply set r0=r1=r, and K0=K1=K.

function y = c_poli_ks_update(P,a_grid,zi,Pr,c_pol,r0,r1,w0,K0,K1)

%Model parameters
beta = P(1);
gamma = P(2);
d = P(4);

%lenghts of the asset share and productivity grids
n_a = length(a_grid);
n_z = length(zi);

c_poli = zeros(n_z,n_a);
%Compute today's consumption using the Euler equation, for fixed values of
%today's savings and productivity.
c_i =  ( ( beta*(1+(r1)-d) ).*( Pr*((1./c_pol).^gamma)  ) ).^(-1/gamma);

%Compute today's implied initial asset shares.
a_grid1 = ( c_i + (a_grid+1)*K1 - (1+zi')*(w0 ) )/(K0*(1+(r0)-d)) - 1;

%Interpolate to complete the mapping from the asset grid to the new
%consumption rule
for j=1:n_z
    
    c_poli(j,:) = interp1(a_grid1(j,:), c_i(j,:), a_grid, 'linear', 'extrap'); 
      
end

%Correct values that violate the constraints (-1<=a_poli<=a_grid(I)):
%Compute the savings rule
a_poli = ( ( 1+zi')*w0 + K0*(1+r0-d)*(a_grid+1) - c_poli)/K1 - 1;
%Consumption if the borrowing constraint binds
c_poli_cd = ( 1+zi')*w0 + K0*(1+r0-d)*(a_grid+1);
%Consumption if the upper limit for savings binds
c_poli_cu = ( 1+zi')*w0 + K0*(1+r0-d)*(a_grid+1) - (a_grid(n_a)+1)*K1;

%Replace values that violate constraints by the proper ones
c_poli(a_poli<=-1)=c_poli_cd(a_poli<=-1);
c_poli(a_poli>=a_grid(n_a))=c_poli_cu(a_poli>=a_grid(n_a));

%Return today's policy rule
y = c_poli;