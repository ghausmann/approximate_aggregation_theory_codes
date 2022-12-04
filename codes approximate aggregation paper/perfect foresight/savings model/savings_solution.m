%This script solves the savings model with heterogeneous agents (HA), one
%bond and perfect foresight, from the paper "Approximate Aggregation Theory
%for Heterogeneous-agent models", by Guillermo Hausmann-Guil. The script
%performs the following tasks:
% - It computes the SS of the model, and provides the agents' optimal rules
% and the SS wealth distribution.
% - It computes the dynamic solution of the model, using the M(1), M(2) and M3)
% approximations.
% - It computes the IRFs of each approximation, and compares them with
% the nonlinear perfect-foresight transition.

clear;

disp('------------------------------------------');
disp('Savings model with perfect foresight');
disp('------------------------------------------');

% Model Parameters
beta = 0.9834; %Impatience
gamma =  2; % CRRA elasticity parameter
phi_ss =1; %SS Borrowing limit
su = 0.25; %Standard deviation of idiosyncratic shocks
rho = 0.75; %Persistence parameter of the borrowing limit process

%Vector of parameters
P = [beta gamma phi_ss rho];

%%
%--------------------------------------------------------------------------
%HA solution: STEADY STATE
%--------------------------------------------------------------------------

%Grid-based algorithm parameters
mytol = 1e-12; %Tolerance for convergence
n_b = 100; %Number of grid points for asset holdings
bmin   = -(phi_ss); % lower bound (make it lower than phi_ss for a positive shock to phi)
bmax   =10; %upper bound
%Asset grid
b_grid = [bmin bmin + ((1:(n_b-1))/(n_b-1)).*(bmax - bmin)];

n_u = 7; %Number of grid points for the idiosyncratic transitory shock
%Use this to approximate the individual iid. shock using the Rouwenhorst's method
[Pr, x_Rouw] = rouwen(0, 0,su,n_u);
ui = x_Rouw'; %Individual income shock grid
Pr = Pr'; %Probability transition matrix (rows add up to one)

%--------------------------------------------------------------------------
%Initial guesses
q_guess = beta*1.01; %SS bond price

b_pol_guess =  0.9*b_grid + 0.9*ui';  %bond rule guess
b_pol_guess(b_pol_guess<=(-phi_ss)) =  b_grid(1);

c_pol = exp(ui') + b_grid - q_guess*b_pol_guess; %consumption rule guess
c_pol(c_pol<=1e-6) = 1e-6;

%Save the guess (so that it can be used by the solver)
save('my_data.mat','c_pol');

%Handle function to solve. It takes as an input a candidate SS bond price,
%computes the optimal consumption rule using the endogenous grid method,
%derives the corresponding stationary wealth distribution, and returns the
%excess demand.
fun_savings_ss = @(x)my_savings_ss(P,mytol,b_grid,ui,Pr,x);
tic
qss = fsolve(fun_savings_ss,q_guess);
toc
rss = 100*(1/qss - 1)

%Compute the optimal consumption rule
dif=1;
while dif > mytol

    %Mapping from the old consumption rule to the new one, using the
    %endogenous grid method.
    c_poli = c_poli_savings_update(P,b_grid,ui,Pr,c_pol,phi_ss,qss);
    dif = max(max(abs(c_poli - c_pol)));  %Check covergence
    c_pol = c_poli; %Replace the old rule with the new one

end

%Compute the corresponding rule for savings
b_poli = (exp(ui') + b_grid - c_pol)/qss;
%Make sure the constraints hold.
b_poli(b_poli<=(-phi_ss) )= -phi_ss;
b_poli(b_poli>=(b_grid(n_b))) = b_grid(n_b);


%Compute the associated stationary distribution (n_u-by-n_b matrix)
D0 = ones(n_u,n_b) / (n_u*n_b); %Initial guess, uniform distribution
Dss = ss_distribution(b_grid,Pr,b_poli,D0,mytol);

mDss = sum(Dss);  %Marginal wealth distribution
muss = sum(Dss,2);  %Marginal income distribution

%Compute moments of the stationary distribution
Bss = sum(sum(Dss.*b_poli)); %Mean
Vss = sum(sum(Dss.*(b_poli.^2))) ;%Variance
M3ss = sum(sum(Dss.*(b_poli.^3))) ;%3rd moment ("skewness")

%%
tic
%--------------------------------------------------------------------------
%HA solution: AGGREGATE DYNAMICS
%--------------------------------------------------------------------------
%Maximum number of periods (by then aggregates must have converged to the SS)
T = 100;
dS = 0.00001;  %Step size for differentiation

%Compute a sequence of average derivatives of functions of today's savings
%rule with respect to the sequences of borrowing limits and bond prices from today until T.
idgdX = savings_idgdX(P,qss,b_grid,ui,Pr,c_pol,Dss,dS,T);
idgdphi = idgdX(1:3,:); %for the borrowing limit
idgdq = idgdX(4:6,:); %for the bond price

%Compute the average partial derivatives of functions of today's
%savings rule with respect to current individual states
idgds = savings_idgds(b_grid,b_poli,Pr,Dss);

%--------------------------------------------------------------------------
%HA local solution, M(1) approximation.
%The recursive solution takes the form: (q-q_ss)  =  lZ*(Z-Z_ss)

%We have a explicit formula for lZ , which uses as inputs the sequences of
%average derivatives.
lZ1 = savings_M1_Z(rho,T,idgdq(1,:),idgdphi(1,:));

%--------------------------------------------------------------------------
% HA local solution, M(2) approximation
% The recursive solution for the equilibrium bond price takes the form
%(q-q_ss) = lV*(V-Vss) + lZ*(Z-Z_ss)
%(Vp-Vss) = mV*(V-Vss) + mZ*(Z-Z_ss)

%We use a system of 2 equations to solve for the coefficients of the
%endogenous states {lV,mV}, using as inputs the average derivatives.
M2_S_fun = @(x)savings_M2_S_fun(T,idgds(1,1:2),idgdq,x);
guess2=[0.5 0.5];
M2_S_sol = fsolve(M2_S_fun,guess2);
%We solve for the coefficients of the exogenous states with basic matrix
%algebra, using as inputs the already solved coefficients, and the
%sequences of average derivatives.
M2_Z_sol = savings_M2_Z(rho,M2_S_sol,T,idgdq,idgdphi);

lV2 = M2_S_sol(1);
mV2 = M2_S_sol(2);

lZ2 = M2_Z_sol(1);
mZ2 = M2_Z_sol(2);

% %--------------------------------------------------------------------------
% HA local solution, M(3) approximation
% The recursive solution for the equilibrium bond price now takes the form
%(q-q_ss)   = lV*(V-Vss) + lM3*(M3-M3ss) + lZ*(Z-Z_ss)
%(Vp-Vss)   = mV*(V-Vss) + mM3*(M3-M3ss) + mZ*(Z-Z_ss)
%(M3p-M3ss) = nV*(V-Vss) + nM3*(M3-M3ss) + nZ*(Z-Z_ss)

%We use a system of 6 equations to solve for the coefficients of the
%endogenous states {lV,lM3,mV,mM3,nV,nM3}.
M3_S_fun = @(x)savings_M3_S_fun(T,idgds,idgdq,x);
guess3 = [lV2 0 mV2 0 0 0];
M3_S_sol = fsolve(M3_S_fun,guess3);

%We solve for the coefficients of the exogenous states with basic matrix
%algebra, using as inputs the already solved coefficients, and the
%sequences of average derivatives.
M3_Z_sol = savings_M3_Z(rho,M3_S_sol,T,idgdq,idgdphi);

lV3 = M3_S_sol(1);
lM3 = M3_S_sol(2);
mV3 = M3_S_sol(3);
mM3 = M3_S_sol(4);
nV3 = M3_S_sol(5);
nM3 = M3_S_sol(6);

lZ3 = M3_Z_sol(1);
mZ3 = M3_Z_sol(2);
nZ3 = M3_Z_sol(3);

toc
%--------------------------------------------------------------------------
%IRFs borrowing limit
%--------------------------------------------------------------------------
%Calculate the IRFs for the M(1), M(2), and M(3) approximations
phi_t = zeros(T,1);

qt1 = zeros(T,1);

qt2 = zeros(T,1);
Vt2 = zeros(T,1);

qt3 = zeros(T,1);
Vt3 = zeros(T,1);
Mt3 = zeros(T,1);

Dz = -0.01; %Shock to the borrowing limit
phi_t(1) = phi_ss + Dz;

qt1(1) = qss + lZ1*Dz;

qt2(1) = qss + lZ2*Dz;
Vt2(1) = Vss+ mZ2*Dz;

qt3(1) = qss + lZ3*Dz;
Vt3(1) = Vss+ mZ3*Dz;
Mt3(1) = M3ss+ nZ3*Dz;

%Compute the IRFs using the recursive formulas
for t=2:T

    phi_t(t) = (1-rho)*phi_ss + rho*phi_t(t-1);
    qt1(t) = qss + lZ1*(phi_t(t)-phi_ss);

    Vt2(t) = Vss + mV2*(Vt2(t-1) - Vss) + mZ2*(phi_t(t)-phi_ss);
    qt2(t) = qss + lV2*(Vt2(t-1) - Vss) + lZ2*(phi_t(t)-phi_ss);


    Vt3(t) = Vss + mV3*(Vt3(t-1) - Vss)    + mM3*(Mt3(t-1) - M3ss)  + mZ3*(phi_t(t)-phi_ss);
    qt3(t) = qss + lV3*(Vt3(t-1) - Vss) + lM3*(Mt3(t-1) - M3ss)  + lZ3*(phi_t(t)-phi_ss);
    Mt3(t) = M3ss + nV3*(Vt3(t-1) - Vss)   + nM3*(Mt3(t-1) - M3ss)  + nZ3*(phi_t(t)-phi_ss);

end

%--------------------------------------------------------------------------
%Accuracy Check
%--------------------------------------------------------------------------

%I check the accuracy of the approximations in two ways. In all cases, I
%use the function "savings_transition.m", which builds the sequence of T
%excess demands for a given sequence of borrowing limits, and a sequence of candidate
%bond prices.
%First, I evaluate the transition function with the three approximations,
%and check the sum of the squared errors (excess demands). The closer to
%zero, the better the approximation.

acc_check1 = savings_transition(P,T,phi_t,b_grid,ui,Pr,c_pol,Dss,qss,qt1);
my_acc_check1 = sum(acc_check1.^2)

acc_check2 = savings_transition(P,T,phi_t,b_grid,ui,Pr,c_pol,Dss,qss,qt2);
my_acc_check2 = sum(acc_check2.^2)

acc_check3 = savings_transition(P,T,phi_t,b_grid,ui,Pr,c_pol,Dss,qss,qt3);
my_acc_check3 = sum(acc_check3.^2)

%Plot the absolute value of the errors over time
%figure;plot(abs(acc_check1),'r');hold on;plot(abs(acc_check2),'g');plot(abs(acc_check3),'k');

%Next, use the same function to solve for the path of equilibrium bond
%prices.
myfullfun = @(x)savings_transition(P,T,phi_t,b_grid,ui,Pr,c_pol,Dss,qss,x);
qt_th = fsolve(myfullfun,qt3); %,optimset('TolFun',1e-8));

mycomp = [qt_th qt1 qt2 qt3];
mynorm1 = norm(qt_th -qt1)
mynorm2 = norm(qt_th -qt2)
mynorm3 = norm(qt_th -qt3)

%Compute the path of the variance of the wealth distribution over time.
Vt_th = variance_transition(P,T,phi_t,b_grid,ui,Pr,c_pol,Dss,qss,qt_th);
%Check accuracy
vcheck2 = norm(Vt_th-Vt2)
vcheck3 = norm(Vt_th-Vt3)

%IRFs for the bond price (% deviation from SS), for each approximation.
pqt_th = 100*(qt_th - qss)/qss;
pqt1 = 100*(qt1 - qss)/qss;
pqt2 = 100*(qt2 - qss)/qss;
pqt3 = 100*(qt3 - qss)/qss;

%IRFs for the variance (% deviation from SS), for each approximation.
pVt_th = 100*(Vt_th - Vss)/Vss;
pVt2 = 100*(Vt2 - Vss)/Vss;
pVt3 = 100*(Vt3 - Vss)/Vss;

Tp = 30;Tpl = 60;
time = 0:Tp;
timel = 0:Tpl;

figure; %Figure 3 in the paper
subplot(1,2,1);
plot(time,pqt_th(1:Tp+1),'b');hold on;plot(time,pqt1(1:Tp+1),'r:');plot(time,pqt2(1:Tp+1),'g--');plot(time,pqt3(1:Tp+1),'k-.');hold off;
subplot(1,2,2);
plot(timel,pVt_th(1:Tpl+1));hold on;plot(timel,pVt2(1:Tpl+1),'g--');plot(timel,pVt3(1:Tpl+1),'k-.');

figure; %Figure 4 in the paper
subplot(1,2,1);plot(b_grid(1:28),b_grid(1:28),'g:');hold on;plot(b_grid(1:28),b_poli(2,1:28),'b');plot(b_grid(1:28),b_poli(4,1:28),'r--');hold off;
subplot(1,2,2);plot(b_grid(1:65),mDss(1:65),'b');