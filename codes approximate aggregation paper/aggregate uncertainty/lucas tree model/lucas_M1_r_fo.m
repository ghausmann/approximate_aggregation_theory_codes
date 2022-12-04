%This function takes as inputs the pair (pss,lY) that characterizes the
%equilibrium asset price as a function of the type: p = pss + lY*(Y-1),
%where Y is output. Using them, it computes the optimal consumption rule
%with the endogenous grid method, derives the corresponding stationary
%distribution, and returns the excess demand (mean of the distribution) and
%the average derivative of the stationary policy rule w.r.t output. We
%reach an approximate equilibrium when both objects equal zero.
%Inputs are P (vector of parameters), mytol (tolerance value for
%convergence), a_grid and zi (grids for asset and productivity values with
%lengths n_a and n_z), Pr (Probability transition matrix with size
%n_z-by-n_z),  Yg (the grid for output of size n_Y), Pr_Y (Probability
%transition matrix of size n_Y-by-n_Y), and x (the coefficients).

function y = lucas_M1_r_fo(P,mytol,a_grid,zi,Pr,Yg,PrY,x)

%Model parameters
alpha = P(3);

%Asset price coefficients
pss = x(1);
lY = x(2)*pss; 

%lenghts of the asset, productivity, and output grids
n_a = length(a_grid);
n_z = length(zi);
n_Y = length(Yg);

%load the matrix summarizing the current candidate for the optimal
%consumption rule (a n_z-by-n_a-by-n_Y array). We will use it as the
%initial guess for the next candidate rule (which depends on (pss,lY) ).

load('my_c_pol_r.mat','c_pol_r');

%Iterate the consumption rule until convergence
dif = 1;
while dif > mytol
    
    %Mapping from the old consumption rule to the new one, using the
    %endogenous grid method.
    c_poli_r = c_poli_lucas_r_update(P,[pss lY],a_grid,zi,Pr,Yg,PrY,c_pol_r,1);
    
    %Check covergence
    dif = max(max(max(abs(c_poli_r - c_pol_r))));
    %Replacle the old rule with the new one
    c_pol_r = c_poli_r;
    
end

%Save the new consumption rule (will be the new guess for the next
%iteration of the solver).
save('my_c_pol_r.mat','c_pol_r');

%Select the steady-state policy function (the case Y = 1 )
sm_Y = median(1:n_Y);
%Compute the corresponding rule for savings
a_poli = ((1+zi')*(1-alpha) + (1+a_grid)*(pss+alpha) - c_pol_r(:,:,sm_Y))/pss - 1;
%Make sure the constraints hold.
a_poli(a_poli<=-1) = -1;
a_poli(a_poli>=(a_grid(n_a))) = a_grid(n_a);

%Compute the associated stationary distribution (another n_z-by-n_a matrix)
D0  = ones(n_z,n_a) / (n_z*n_a);
Dss = ss_distribution(a_grid,Pr,a_poli,D0,mytol);

%Compute the mean of the distribution
Ap = sum(sum(Dss) .* a_grid);

%Compute the average derivative of the stationary asset rule w.r.t. output.
dAp_dY = lucas_M1_r_dgdY_fo(P,a_grid,zi,Yg-1,c_pol_r,Dss,0.001,[pss lY]);

%Return the equations
y = [Ap;dAp_dY];
 
    