%This code computes the first-order solution of the Lucas Tree model with
%aggregate risk, from the paper "Approximate Aggregation Theory for
%Heterogeneous-agent models", by Guillermo Hausmann-Guil.
%The script performs the following tasks:
%- It computes the perfect foresight solution
%(steady state + aggregate dynamics)
%- It computes the first-order solution with aggregate
%- It runs the den Haan (2010) accuracy test


clear;

disp('-------------------------------------');
disp('Lucas tree model with aggregate risk');
disp('-------------------------------------');

%Model parameters
beta = 0.99; %Impatience
gamma = 2; % CRRA elasticity parameter
%gamma = 1 %for exact solution
alpha = 0.1; % Share of tree dividends from total GDP
sY = 0.01; %std of output shocks
%sY = 0.02; %std for High risk scenario
rho = 0.95; %persistence of output process
%Unconditional Standard deviation
std_Y = sY/((1-rho^2)^0.5);

%Vector of parameters
P = [beta gamma alpha sY rho];

%--------------------------------------------------------------------------
%RA solution
%--------------------------------------------------------------------------

%RA local solution, for comparison purposes and initial guess
%It takes the form p = pss_RA + ly1_RA *(Y-1)
%RA Asset price at SS
pss_RA = alpha*beta/(1-beta);
%RA Marginal response to Y (first and second order)
ly1_RA = -( beta*gamma*(alpha + pss_RA)*(rho - 1) - alpha*beta*rho )/(1 - beta*rho);
%In ratios
l1_RA = ly1_RA/pss_RA;

%%
%--------------------------------------------------------------------------
%HA solution: STEADY STATE
%--------------------------------------------------------------------------

%Grid-based algorithm parameters
mytol = 1e-12; %Tolerance for convergence
n_a = 100; %Number of grid points for asset holdings
amin   = -1;  %lower bound
amax   = 4; %upper bound
%Asset grid
a_grid = [amin amin + ((1:(n_a-1))/(n_a-1)) * (amax - amin)];

%--------------------------------------------------------------------------
%Idiosincratic process (same as in KS model)
zi = [0 1] %Individual labor states

Peu = 0.025; %job separation rate
pi_ss = 0.95; %SS employment target
Pue = (pi_ss/(1-pi_ss))*Peu; % job finding rate
Pr = [(1-Pue) Pue;Peu (1-Peu)] %Probability transition matrix (rows add up to one)
n_z = length(zi); %number of gridpoints for labor

%Iterate over the probability transition matrix to compute the
%unconditional distribution, and use the mean to normalize the process
dif=1;
mp0 = ones(n_z,1)/n_z;
while dif>mytol
    
    mp1 = Pr'*mp0;
    dif = norm(mp1-mp0);
    mp0 = mp1;
    
end
mmean0 = sum(zi.*mp0'); %compute the unconditional mean
zi = zi/mmean0 - 1; %normalize the idosincratic process
m_xi = sum(zi.*mp0'); %check normalization

%--------------------------------------------------------------------------
%Output process
n_Y = 5; %Number of grid points for aggregate output (make sure it is odd!)
%Use this to approximate the AR(1) process using the Rouwenhorst's method
[PrY, Z_Rouwz] = rouwen(rho, 0,std_Y,n_Y);
Yg = exp(Z_Rouwz)'; %output grid
PrY = PrY'; %Probability transition matrix (rows add up to one)

%--------------------------------------------------------------------------
%Initial guesses
p_guess = pss_RA*1.1; %asset price guess
a_pol_guess =  0.9*(a_grid) + zi'; %asset rule guess
a_pol_guess(a_pol_guess<=-1) =  -1;

c_pol = (1+zi')*(1-alpha) + (pss_RA+alpha)*(1+a_grid) - pss_RA*(1+a_pol_guess); %consumption rule guess
c_pol(c_pol<=1e-6) = 1e-6;

%Save the guess (so that it can be used by the solver)
save('my_data.mat','c_pol');

%Handle function to solve. It calls the function "my_lucas_ss.m", which
%takes as an input a candidate SS asset price, computes the optimal
%consumption rule using the endogenous grid method, derives the
%corresponding stationary wealth distribution, and returns the excess demand
%(mean of the distribution minus 1). See the function for details.
fun_lucas_ss = @(x)my_lucas_ss(P,mytol,a_grid,zi,Pr,x);
tic
%Finds the SS equilibrium price using the Matlab non-linear solver
pss_d = fsolve(fun_lucas_ss,p_guess);
toc

%Load back the optimal consumption rule
load('my_data.mat','c_pol');

%Re-compute the optimal consumption rule (for accuracy purposes)
dif=1;
while dif > mytol
    
    %Mapping from the old consumption rule to the new one, using the
    %endogenous grid method (see the function "c_poli_lucas_update.m" for
    %details).
    c_poli = c_poli_lucas_update(P,a_grid,zi,Pr,c_pol,1,1,pss_d,pss_d);
    dif = max(max(abs(c_poli - c_pol)));
    c_pol = c_poli;
    
end

%Compute the corresponding rule for savings
a_poli = ( (1+zi')*(1-alpha) + (1+a_grid)*(pss_d+alpha) - c_pol)/pss_d - 1;
%Make sure the constraints hold.
a_poli(a_poli<=-1) = -1;
a_poli(a_poli>=(a_grid(n_a))) = a_grid(n_a);

%Compute the stationary wealth distribution
D0 = ones(n_z,n_a) / (n_z*n_a); %Initial guess, uniform distribution
Dss = ss_distribution(a_grid,Pr,a_poli,D0,mytol);

%Compute and plot the marginal stationary distribution of assets
mD = sum(Dss);
figure;plot(a_grid,mD);
% 
%%
%--------------------------------------------------------------------------
%AGGREGATE DYNAMICS: perfect foresight
%--------------------------------------------------------------------------
%Maximum number of periods (by then aggregates must have converged to the SS)
T = 100;
%Step size for differentiation
dS = 0.00001;

tic
%Compute a sequence of average derivatives of today's asset rule with
%respect to the sequences of output and the asset price, from today until
%T.
idgdX = lucas_idgdX(P,pss_d,a_grid,zi,Pr,c_pol,Dss,dS,T);

%HA local solution, M(1) approximation.
%As in the RA case, the recursive solution takes the form
%p = p_ss + lY*(Y-1)
%where p_ss is the SS asset price, Y aggregate output, and ly
%the marginal response of the asset price to changes in Y. In this case we
%have a explicit theoretical formula, which uses as inputs the sequences of
%average derivatives.
lYd = lucas_M1_pf(rho,T,idgdX);

toc
%%
%--------------------------------------------------------------------------
%AGGREGATE DYNAMICS: aggregate risk
%--------------------------------------------------------------------------

%Initial guess for the consumption rule (an SxIxSz array)
c_pol_r = c_pol.*ones(n_z,n_a,n_Y);
%Save the guess (so that it can be used by the solver)
save('my_c_pol_r.mat','c_pol_r');

%FIRST-ORDER SOLUTION: p = pss_r + lYr*(Y-1)
%Handle function to solve for the coefficients. 
%Note: the marginal response enters as a ratio wrt. SS price
myfun_risk = @(x)lucas_M1_r_fo(P,mytol,a_grid,zi,Pr,Yg,PrY,x);  
%Solve for the approximate equilibrium using the Matlab non-linear solver
tic
sol_risk = fsolve(myfun_risk,[pss_d l1_RA]) 
toc

pss_r = sol_risk(1);
lYr= sol_risk(2)*pss_r;

load('my_c_pol_r.mat','c_pol_r');

%Iterate the consumption rule until convergence
dif = 1;
while dif > mytol
    
    %Mapping from the old consumption rule to the new one, using the
    %endogenous grid method (see the function "c_poli_lucas_r_update.m" for
    %details).
    c_poli_r1 =  c_poli_lucas_r_update(P,[pss_r lYr],a_grid,zi,Pr,Yg,PrY,c_pol_r,1);
                   
    %Check covergence
    dif = max(max(max(abs(c_poli_r1 - c_pol_r))));
    %Replacle the old rules with the new ones
    c_pol_r = c_poli_r1;
    
end
%Select the stationary policy function (associated with z = (Y-1) = 0 )
sm_Y = median(1:n_Y);
%Compute the corresponding rules for savings
a_poli_r = ((1+zi')*(1-alpha) + (1+a_grid)*(pss_r+alpha) - c_pol_r(:,:,sm_Y))/pss_r - 1; %linear case

%Make sure the constraints hold.
a_poli_r(a_poli_r<=-1) = -1;
a_poli_r(a_poli_r>=(a_grid(n_a))) = a_grid(n_a);

%Compute the corresponding stationary wealth distributions
Dss_r = ss_distribution(a_grid,Pr,a_poli_r,D0,mytol);

%%
%--------------------------------------------------------------------------
%Accuracy test
%--------------------------------------------------------------------------

%First, perform a basic check, evaluating the function "my_check_r.m",
%which returns the excess demand for a given current price (here the steady
%state price pss_r), and a given level of output (here Y=1), using the
%approximate solution to form expectation about future prices.
basic_check =  my_check_r(P,a_grid,zi,Pr,Yg,PrY,c_pol_r,[pss_r lYr],Dss_r,sm_Y,pss_r,1)

%Next, perform the den Haan (2010) test, by computing absolute percent
%errors.
Tsim = 2000; %length of the simulated time series
Dt = Dss_r;  %current wealth distribution
pt = zeros(1,Tsim); %vector for current asset price
pa = zeros(1,Tsim); %vector for approximate asset price
myflag = zeros(1,Tsim); 

myrand = rand(Tsim,1);
%generate a time series using the Rouwenhorst method
my_ind = hitm_s(PrY',myrand); %indices for Yg
%uncomment to get the numbers from Table 4
%load('my_hit_baseline.mat','my_ind');
%load('my_hit_highrisk.mat','my_ind');

%Perform the test
for t = 1:Tsim
    
    ind_Y = my_ind(t); %pick the current level of output
    %compute the approximate asset price
    pa(t) = pss_r + lYr*(Yg(ind_Y)-1); 
     
    %Handle functions to solve
    fun_check1 = @(x)my_check_r(P,a_grid,zi,Pr,Yg,PrY,c_pol_r,[pss_r lYr],Dt,ind_Y,x,1);
    %Solve for the current asset price that clears the market
    [pt(t),~,myflag(t)]  = fzero(fun_check1,pa(t));
    
    %Update the wealth distribution, using the function
    %"my_output_fo.m".
    [~,Dt] = my_output_r(P,a_grid,zi,Pr,Yg,PrY,c_pol_r,[pss_r lYr],Dt,ind_Y,pt(t),1);
            
end

pt = pt(1001:Tsim);
pa = pa(1001:Tsim);

%Compute the errors, and return statistics and histogram
pct_errors1=abs(100*(pt-pa)./pa); 
mean_errors1 = mean(pct_errors1); max_errors1 = max(pct_errors1);
Statistics = [mean_errors1 max_errors1]
figure;histogram(pct_errors1); title('linear');
