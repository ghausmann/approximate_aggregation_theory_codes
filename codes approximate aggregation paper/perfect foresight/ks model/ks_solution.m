%This script solves a version of the Krusell-Smith model with perfect
%foresight, from the paper "Approximate Aggregation Theory for
%Heterogeneous-agent models", by Guillermo Hausmann-Guil. The script
%performs the following tasks:
% - It computes the steady-state (SS) of the model, and provides the
% agents' optimal rule and the SS wealth distribution.
% - It computes the dynamic solution to the model, using the M(1)
% approximation.
% - It computes the IRF of the equilibrium rental price, and compares it
% with the "theoretical" IRF (computed using the Implicit function theorem
% around the SS).

clear;

disp('------------------------------------------');
disp('Krusell-Smith model with perfect foresight');
disp('------------------------------------------');


% Model Parameters
beta = 0.99; %Impatience
gamma =1; % CRRA elasticity parameter
alpha = 0.3; % Cobb-Douglas elasticity parameter
d = 0.025; % Depreciation rate
rho = 0.95; %Persistence parameter of the aggregate productivity process

%Vector of model parameters
P = [beta gamma alpha d rho];

%Grid-based algorithm parameters
mytol = 1e-12; %Tolerance for convergence
n_a = 200; %Number of grid points for shares of asset holdings
amin   = -1;  %lower bound
amax   = 3; %upper bound
T = 300; %Maximum number of periods (by then aggregates must have converged to the SS)

%SS from RA model: aggregate capital and associated prices
Kss_ra= ((1/beta - (1-d))/alpha)^(-1/(1-alpha));
wss_ra = ((1-alpha))*(1-alpha)*Kss_ra^alpha;
rss_ra = alpha*(Kss_ra)^(alpha-1);
Yss_ra = ((1-alpha))*Kss_ra^alpha;

%%
%--------------------------------------------------------------------------
%HA solution: STEADY STATE (SS)
%--------------------------------------------------------------------------

%Asset share grid
a_grid = [amin amin + ((1:(n_a-1))/(n_a-1)) * (amax - amin)];

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
m_zi = sum(zi.*mp0'); %check normalization

%--------------------------------------------------------------------------
%Initial guesses
r_guess = 0.99*rss_ra; %SS rental rpice guess
K_guess = (alpha/r_guess)^(1/(1-alpha));
a_pol_guess =  0.99*a_grid + 0.6*(zi)'; %asset rule guess
a_pol_guess(a_pol_guess<=-1) =  0;
c_pol = (1+zi')*wss_ra + K_guess*(1+rss_ra-d)*a_grid - K_guess*a_pol_guess; %consumption rule guess
c_pol(c_pol<=1e-6) = 1e-6;

%Save the guess (so that it can be used by the solver)
save('my_data.mat','c_pol');

%Handle function to solve. It  takes as an input a candidate SS rental
%price, computes the optimal consumption rule using the endogenous grid
%method, derives the implied stationary wealth distribution, and returns
%the excess demand. 
fun_ks_ss = @(x)my_ks_ss(P,mytol,a_grid,zi,Pr,x);
tic
%Find the SS equilibrium rental price using the Matlab non-linear solver
rss_ha = fsolve(fun_ks_ss,r_guess,optimset('TolFun',1e-5))
toc

%Implied SS aggregate capital, wage rate, and interest rate
Kss_ha = ((alpha/rss_ha)^(1/(1-alpha)));
wss_ha = (1-alpha)*(Kss_ha^alpha)*((1-alpha));
irate_ha = 100*(rss_ha-d)

%Partial derivatives of the rental price w.r.t. today's aggregate capital
%and productivity, evaluated at the SS. They will be useful later for
%computing IRFs.
rK = -(1-alpha)*rss_ha/Kss_ha;
rZ =rss_ha;

%Compute the optimal consumption rule, iterating until convergence
dif=1;
while dif > mytol

    %Compute the new consumption rule
    c_poli = c_poli_ks_update(P,a_grid,zi,Pr,c_pol,rss_ha,rss_ha,wss_ha,Kss_ha,Kss_ha);
    dif = max(max(abs(c_poli - c_pol)));  %Check covergence
    c_pol = c_poli; %Replace the old rule with the new one

end

%Compute the corresponding rule for savings
a_poli = ((1+zi)'*wss_ha + Kss_ha*(1+rss_ha-d)*(a_grid+1) - c_pol)/Kss_ha - 1;

%Make sure the constraints hold.
a_poli(a_poli<=-1) = -1;
a_poli(a_poli>=(a_grid(n_a))) = a_grid(n_a);

%Compute the associated stationary distribution Dss (n_a-by_n_z matrix)
D0 = ones(n_z,n_a) / (n_z*n_a); %Initial guess, uniform distribution
Dss = ss_distribution(a_grid,Pr,a_poli,D0,mytol);
mDss = sum(Dss); %Marginal wealth distribution
%Compute mean of the stationary distribution and check a zero
%excess-demand.
Mss1 = sum(mDss.*(a_grid));
check_ss = Mss1

%%
%--------------------------------------------------------------------------
%HA solution: AGGREGATE DYNAMICS using the M(1) approximation
%--------------------------------------------------------------------------

% The M(1) approximation treats aggregate capital as an endogenous
% aggregate state, with law of motion
%Kp = K_ss + lK*(K-Kss) + lZ*(Z-1)
%where K is the current level of aggregate capital (Kss its SS value), Z
%current aggregate productivity, with exogenous process
%Zp = 1 + rho*(Z-1)
%and lK and lZ two coefficients to be determined,  for which I solve next.

tic

Tmax = T;
dS = 0.00001; %step size for differentiation

%Compute a sequence of average derivatives of functions of today's savings rule with respect to the sequence
%of prices and capital, from today until T.
idgdX = ks_idgdX(P,rss_ha,a_grid,zi,Pr,c_poli,Dss,dS,Tmax);

% %Uncomment to compute the coefficient lK using the Matlab non-linear solver
% lKg = 0.5;
% funlK = @(x)ks_M1_S_fun(P,rss_ha,T,idgdX,x);
% lK = fsolve(funlK,lKg);

%Find lK using time iteration
lKg = unifrnd(-0.99,0.99); %random initial guess
mdiff = 1;
s1 = 1;
cM10 = lKg;

while mdiff>1e-9

    %See the function ks_M1_S_fun_it.m for details
    cM11 = ks_M1_S_fun_it(P,rss_ha,Tmax,idgdX,dS,cM10);
    mdiff = (cM11-cM10)^2;
    cM10 = cM11;
    s1 = s1+1;

end
lK=cM11;

%Compute lZ using the formula from the paper.
lZ = ks_M1_Z(P,rss_ha,T,idgdX,lK);

toc

M1_solution = [lK lZ]

%%
%--------------------------------------------------------------------------
%IRFs
%--------------------------------------------------------------------------
%Calculate the IRFs for aggregate capital and rental price.

Zt = ones(1,T);
Kt = Kss_ha*ones(1,T);
rt = rss_ha*ones(1,T);

dZ0 = 0.00001; %small shock to aggregate productivity
Zt(1) = 1 + dZ0;
rt(1) = rss_ha + rZ*(Zt(1)-1);
Kt(1) = Kss_ha + lZ*(Zt(1)-1);

for t=2:T

    Zt(t) = 1 + rho*(Zt(t-1)-1);
    Kt(t) = Kss_ha + lK*(Kt(t-1)-Kss_ha) + lZ*(Zt(t)-1);
    rt(t) = rss_ha + rK*(Kt(t-1)-Kss_ha) + rZ*(Zt(t)-1);

end

%%
%--------------------------------------------------------------------------
%Accuracy Check
%--------------------------------------------------------------------------

%I check the accuracy of the approximation in two ways. In both cases, I
%use the function "my_ks_transition.m", which computes a sequence of T
%excess demands for a given productivity shock. See the function for
%details.

%First, I evaluate the transtition function with my M(1) approximation, and
%check the sum of the squared errors. The closer to zero, the better the
%approximation.
rt_app = rt(2:(T));
acc_check = my_ks_transition(P,T,dZ0,a_grid,zi,Pr,c_pol,rss_ha,Dss,rt_app);
my_check = sum(acc_check.^2)
figure;plot(acc_check );

%Second, I compute the Jacobian of the system of equations formed by the
%sequence of T market-clearing conditions, and use it to compute the
%theoretical IRF (theoretical derivatives of rental prices w.r.t. a current
%productivity shock) following the Implicit Function theorem. Then I
%compare it with my approximation.
tic
vrss=rss_ha*ones(1,T-1);
%Calculate the Jacobian, using the fsolve command ("Direct method" in Auclert
%et al. (2021))
funtrans_ss = @(x)my_ks_transition(P,T,0,a_grid,zi,Pr,c_pol,rss_ha,Dss,x);
[~,F0,~,~,Jt] = fsolve(funtrans_ss,vrss,optimset('Tolfun',1e-5));
%Evaluate the transition function at the SS prices with the shock. We use
%this to construct the partial derivatives of the sequence of functions
%w.r.t. the shock:
F1 = my_ks_transition(P,T,dZ0,a_grid,zi,Pr,c_pol,rss_ha,Dss,vrss);
dZt = -(F1-F0)/dZ0;
toc

% save('my_results.mat','rss_ha','Jt','F0','F1','dZt');

%Theoretical IRF and sequence of interest rate following a 1% shock to
%productivity Z0.
irf_th = (1/rss_ha)*[rZ ((Jt^-1)*dZt')'];
rt_th = 100*(rss_ha + irf_th*(rss_ha/100) - d);

%Same for the from M(1) approximation.
irf_M1 = (1/rss_ha)*(rt-rss_ha)/dZ0;
rt_M1 = 100*(rss_ha + irf_M1*(rss_ha/100) - d);

mycomp = [irf_th ;irf_M1];
mynorm1 = norm(irf_th - irf_M1) %Distance check
diff_M1 = 100*abs(rt_M1-rt_th)./rt_th; %percent deviation check

%Plot the results
Tp = T;
time = 0:Tp-1;
figure; %Figure 1 in the paper
subplot(1,2,1);plot(time,rt_th(1:Tp ),'b');hold on;plot(time,rt_M1(1:Tp  ),'r--');hold off;
subplot(1,2,2);plot(time,diff_M1(1:Tp  ),'r--');

figure; %Figure 2 in the paper
subplot(1,2,1);plot(a_grid(1:51),a_grid(1:51),'g:');hold on;plot(a_grid(1:51),a_poli(1,1:51),'b');plot(a_grid(1:51),a_poli(2,1:51),'r--');hold off;
subplot(1,2,2);plot(a_grid,mDss,'b');
