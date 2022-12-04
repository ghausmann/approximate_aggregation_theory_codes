%This function takes as an input a candidate SS rental price of capital,
%computes the optimal consumption rule using the endogenous grid method,
%derives the corresponding invariant distribution, and returns the excess
%demand.
% Inputs are P (vector of parameters), mytol (tolerance value for
%convergence), a_grid and xi (grids for asset shares and productivity
%values with lengths n_a and n_z), Pr (Probability transition matrix with
%size n_z-by-n_z), and r0 (the candidate rental price).

function y = my_ks_ss(P,mytol,a_grid,zi,Pr,r0)

%Model parameters
alpha = P(3);
d = P(4);

%lenghts of the asset and productivity grids
n_a = length(a_grid);
n_z = length(zi);

%load the matrix summarizing the current candidate for the optimal
%consumption rule (a n_z-by-n_a matrix). We will use it as the initial guess for
%the next candidate matrix (which depend on the candidate SS rental price
%of capital)
load('my_data.mat','c_pol');

%Implied SS aggregate capital and wage rate
K0 = (alpha/r0)^(1/(1-alpha));
w0 = ((1-alpha))*(1-alpha)*K0^alpha;

%Iterate the consumption rule (a n_z-by-n_a matrix) until convergence
dif = 1;
while dif>mytol

    %Compute the new consumption rule with the function
    %c_poli_ks_update.m (it uses the endogenous gridpoints method)
    c_poli = c_poli_ks_update(P,a_grid,zi,Pr,c_pol,r0,r0,w0,K0,K0);
    dif = max(max(abs(c_poli - c_pol))); %Check covergence
    c_pol = c_poli; %Replace the old rule with the new one

end

%Construct the implied asset share rule
a_poli = ((1+zi')*w0 + K0*(1+r0-d)*(a_grid+1) - c_pol)/K0 - 1;
%Make sure the constraints hold
a_poli(a_poli<=-1) = -1;
a_poli(a_poli>=(a_grid(n_a))) = a_grid(n_a);

%Call the function "ss_distribution.m" to compute the associated stationary
%distribution (another n_z-by-n_a matrix)
D0 = ones(n_z,n_a) / (n_z*n_a); %Initial guess (uniform distribution)
Dss = ss_distribution(a_grid,Pr,a_poli,D0,mytol);
%Compute the mean of the distribution (aggregate demand)
Ys = sum(sum(Dss.*a_poli));

%Save the new consumption rule (will be the new guess for the next
%iteration of the solver).
save('my_data.mat','c_pol');

y = Ys;