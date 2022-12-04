%This function takes as an input a candidate SS asset price, computes the
%optimal consumption rule using the endogenous grid method, derives the
%corresponding invariant distribution, and returns the excess demand (mean
%of the distribution minus 1). The SS equilibrium price will be the one
%which makes the excess demand equal to zero.
%Inputs are P (vector of parameters), mytol (tolerance value for
%convergence), a_grid and xi (grids for asset and productivity values with
%lengths n_a and n_z), Pr (Probability transition matrix with size n_z-by-n_z), and p
%(the asset price).

function y = my_lucas_ss(P,mytol,a_grid,zi,Pr,p)

%Model parameters
alpha = P(3);

%lenghts of the asset and productivity grids
n_a = length(a_grid);
n_z = length(zi);

%load the matrix summarizing the current candidate for the optimal
%consumption rule (an SxI matrix). We will use it as the initial guess for
%the next candidate matrix (which depend on the candidate SS asset price
%p).
load('my_data.mat','c_pol');

%Iterate the consumption rule until convergence
dif = 1;
while dif > mytol
    
    %Mapping from the old consumption rule to the new one, using the
    %endogenous grid method (see the function "c_poli_lucas_update.m" for
    %details).
    c_poli = c_poli_lucas_update(P,a_grid,zi,Pr,c_pol,1,1,p,p);
    
    %Check covergence
    dif = max(max(abs(c_poli - c_pol)));
    %Replacle the old rule with the new one
    c_pol = c_poli;
    
end

%Compute the corresponding rule for savings
a_poli = ( (1+zi')*(1-alpha) + (1+a_grid)*(p+alpha) - c_poli)/p - 1;
%Make sure the constraints hold.
a_poli(a_poli<=-1) = -1;
a_poli(a_poli>=(a_grid(n_a))) = a_grid(n_a);

%Call the function "ss_distribution.m" to compute the associated stationary
%distribution (another SxI matrix)
D0 = ones(n_z,n_a) / (n_z*n_a); %Initial guess (uniform distribution)
Dss = ss_distribution(a_grid,Pr,a_poli,D0,mytol);
%Compute the mean of the distribution
Ap = sum(sum(Dss) .* a_grid);

%Save the new consumption rule (will be the new guess for the next
%iteration of the solver).
save('my_c_pol.mat','c_pol');

%Return the asset excess demand
y = Ap;
