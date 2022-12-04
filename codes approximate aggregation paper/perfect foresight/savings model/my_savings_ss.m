%This function takes as an input a candidate SS bond price q, computes the
%optimal consumption rule using the endogenous grid method, derives the
%corresponding invariant distribution, and returns the excess demand (mean
%of the distribution). The SS equilibrium price will be the one which makes
%the excess demand equal to zero.
%Inputs are P (vector of parameters), mytol (tolerance value for
%convergence), b_grid and ui (grids for asset and income shock values with
%lengths n_b and n_u), Pr (Probability transition matrix with size
%n_u-by-n_u), and q.

function y = my_savings_ss(P,mytol,b_grid,ui,Pr,q)

%Model parameters
phi_ss = P(3);

%lenghts of the asset and productivity grids
n_b = length(b_grid);
n_u = length(ui);

%load the matrix summarizing the current candidate for the optimal
%consumption rule (a n_u-by-n_b matrix). We will use it as the initial guess for
%the next candidate matrix (which depend on the candidate SS bond price q).
load('my_data.mat','c_pol');

%Iterate the consumption rule (a SxI matrix) until convergence
dif = 1;
while dif > mytol

    %Mapping from the old consumption rule to the new one, using the
    %endogenous grid method. 
    c_poli = c_poli_savings_update(P,b_grid,ui,Pr,c_pol,phi_ss,q);
    dif = max(max(abs(c_poli - c_pol))); %Check covergence
    c_pol = c_poli; %Replacle the old rule with the new one

end

%Compute the corresponding rule for savings
b_poli = (exp(ui') + b_grid - c_pol)/q;
%Make sure the constraints hold.
b_poli(b_poli<=(-phi_ss)) = -phi_ss;
b_poli(b_poli>=(b_grid(n_b))) = b_grid(n_b);

%Call the function "ss_distribution.m" to compute the associated stationary
%distribution (another n_u-by-n_b matrix)
D0 = ones(n_u,n_b) / (n_u*n_b); %Initial guess (uniform distribution)
Dss = ss_distribution(b_grid,Pr,b_poli,D0,mytol);
%Compute the mean
Bp = sum(sum(Dss.*b_poli));
%Save the new consumption rule (will be the new guess for the next
%iteration of the solver).
save('my_data.mat','c_pol');

%Return the asset excess demand
y = Bp;

