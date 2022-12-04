%This code solves the Two-asset model with aggregate risk from the paper
%"Approximate Aggregation Theory for Heterogeneous-agent models", by
%Guillermo Hausmann-Guil.
%The script performs the following tasks:
%- It computes the local solution
%- It computes Euler equation errors for different values of the exogenous
%states
%- It computes IRFs following an idiosyncratic innovation

clear;

disp('--------------------------------------');
disp('Two-asset model with aggregate risk');
disp('--------------------------------------');

%Model parameters
beta = 0.99; %Impatience
gamma = 1; % CRRA elasticity parameter
phi = 0.25;% borrowing limit
su = 0.1; %std of individual shocks
sz = 0.01; %std of aggregate shocks
rz = 0.9; %persistence of aggregate productivity
std_z = sz/((1-rz^2)^0.5); %unconditional Standard deviations

%Vector of parameters
P = [beta gamma phi su sz rz];


%%
%--------------------------------------------------------------------------
%Local solution
%--------------------------------------------------------------------------

%Compute the simulated time series for the idiosyncratic innovations
% T = 11000;
% xt = su*normrnd(0,1,[1 T]);
% save('myxt.mat','xt');
load('myxt.mat','xt');
T = length(xt);

%Initial guess
guess = [beta*1.01 0.5 beta*1.01 -0.5 0 0 -0.1 0 0 0.5 0 0];

%solve the system of equations derived in the paper
%see local_2assets_system.m for details.
myfun = @(v)local_2assets_system(P,v,xt);
tic
mysol = fsolve(myfun,guess)
toc

%Solution coefficients
qss = mysol(1);
lz = mysol(2);
pss = mysol(3);
nz = mysol(4);
gb = mysol(5);
ga = mysol(6);
gu = mysol(7);
hb = mysol(8);
ha = mysol(9);
hu = mysol(10);
gss = mysol(11);
hss = mysol(12);

%Check the eigenvalues (stability)
Meig = [ga gb;ha hb];
eig(Meig)
myeig = max(abs(eig(Meig)))

%Check Euler-equation errors when evaluated at the steady-state
myerrors = log10(local_2assets_errors(P,mysol,[gss hss 0 0]))

%%
%--------------------------------------------------------------------------
%Euler equation errors
%--------------------------------------------------------------------------
%
Ie = 201;
umine   = -2*su;  % lower bound
umaxe   = 2*su;   % upper bound

zmine   = -2*std_z;  % lower bound
zmaxe   = 2*std_z;   % upper bound

%Grids for exogenous states ut and zt
due = (umaxe-umine)/(Ie-1);
dze = (zmaxe-zmine)/(Ie-1);
u_gride = umine:due:umaxe;
z_gride = zmine:dze:zmaxe;

errors_u  =zeros(2,Ie);
errors_z  =zeros(2,Ie);

%Compute the Euler-equation errors for different grid values
for t = 1:Ie
    
    errors_u(:,t) = log10(local_2assets_errors(P,mysol,[gss hss u_gride(t) 0]));
    errors_z(:,t) = log10(local_2assets_errors(P,mysol,[gss hss 0 z_gride(t)]));
    
end

%Plot the results
figure;
subplot(1,2,1);plot(u_gride,errors_u(2,:),'b');title('(a) for idiosyncratic shocks');
subplot(1,2,2);plot(z_gride,errors_z(2,:),'r');title('(b) for aggregate shocks');


%%
%--------------------------------------------------------------------------
%IRFs ui shock
%--------------------------------------------------------------------------

S = 101;

btt = zeros(1,S); %safe asset
att = zeros(1,S); %risky asset
ct = zeros(1,S); %consumption

%Initial values following a innovation
att(1) = gss +su*gu;
btt(1) = hss +su*hu;
ct(1) = 1 + su  + hss + gss - ( pss*att(1)  ) - ( qss*btt(1) ) ;

for s=2:S
    
    att(s) = gss + gb*(btt(s-1) - hss) + ga*(att(s-1) - gss);
    btt(s) = hss + hb*(btt(s-1) - hss) + ha*(att(s-1) - gss);
    ct(s) = 1  + att(s-1) + btt(s-1) - ( pss*att(s)  ) - ( qss*btt(s) ) ;
    
end

%Individual wealth
wtt = qss*btt+pss*att;
time = 0:(S-1);

%Plot the IRFs
figure;
subplot(2,2,1);
plot(time,btt);title('(a) Bond');
axis([time(1) time(S-1) min(btt)*1.2 max(btt)*1.01]);
subplot(2,2,2);
plot(time,att);title('(b) Risky asset');
axis([time(1) time(S-1) min(att)*0.99 max(att)*1.15]);
subplot(2,2,3);
plot(time,wtt);title('(c) Wealth');
axis([time(1) time(S-1) min(wtt)*1.2 max(wtt)*1.01]);
subplot(2,2,4);
plot(time,ct);title('(d) Consumption');
axis([time(1) time(S-1) min(ct)*0.999 max(ct)*1.001]);
