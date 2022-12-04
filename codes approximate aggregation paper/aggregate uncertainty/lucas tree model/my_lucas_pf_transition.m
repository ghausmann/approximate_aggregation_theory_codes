%This function computes a sequence of size T of excess demands for each
%period asset market, given a sequence of candidate asset prices and the
%initial excess output. The equilibrium sequence of prices is the one that
%clears all markets (makes all the T excess demands 0).
%The inputs are P (vector of model parameters), T (number of periods), dD
%(initial shock), the grids a_grid (with length n_a) and zi (with length
%n_z), Pr (Probability transition matrix with size n_z-by-n_z), the SS
%consumption rule (a n_z-by-n_a matrix),the SS asset price pss, the
%stationary wealth distribution Dss (another n_z-by-n_a matrix), and the
%sequence of asset prices x.

function y  = my_lucas_pf_transition(P,T,dS,a_grid,zi,Pr,c_pol,pss,Dss,x)

%Model parameters
alpha = P(3);
rho = P(5);

n_z = length(zi);
n_a = length(a_grid);

%Sequence of prices + the terminal SS value
pt = [x pss];
%Sequence of excess output
ZZ = zeros(1,T+1);
ZZ(1) = dS;

%Compute the sequence of excess output
for t=2:T
    
    ZZ(t) = rho*ZZ(t-1);
        
end

%--------------------------------------------------------------------------
%BACKWARD ITERATION (of policy rules)
%--------------------------------------------------------------------------

%SS asset rule
a_poli = ( (1+zi')*(1-alpha) + (1+a_grid)*(pss+alpha) - c_pol)/pss - 1;
%Make sure the constraints hold.
a_poli(a_poli<=-1) = -1;
a_poli(a_poli>=(a_grid(n_a))) = a_grid(n_a);

%Compute the sequence of policy functions from the Terminal period T to the
%current period t0, using backward iteration
c_pol_dD = c_pol;

a_polaa = zeros(n_z,n_a,T+1);
a_polaa(:,:,T+1) = a_poli;


for n=(T):-1:1
    
    %Use the function "c_poli_lucas_update.m" to iterate backwards (based
    %on the endogenous grid method)
    c_poli_dD = c_poli_lucas_update(P,a_grid,zi,Pr,c_pol_dD,1+ZZ(n),1+ZZ(n+1),pt(n),pt(n+1));
    %Compute the correspoding savings rule
    a_poli_dD = ( (1+zi')*(1-alpha)*(1+ZZ(n)) + (1+a_grid)*((pt(n) )+alpha*(1+ZZ(n))) - c_poli_dD)/(pt(n) ) -1 ;
    a_poli_dD(a_poli_dD<=-1) = -1;
    a_poli_dD(a_poli_dD>=(a_grid(n_a))) = a_grid(n_a);
    
    c_pol_dD = c_poli_dD;
    a_polaa(:,:,n) = a_poli_dD;
end


%--------------------------------------------------------------------------
%FORWARD ITERATION (of wealth distribution)
%--------------------------------------------------------------------------

%Sequence of wealth distributions
Dt = zeros(n_z,n_a,T+1);
%Sequence of market-clearing conditions
market_cleart = zeros(T,1);

Dt(:,:,1) = Dss;

%Iterate forward
for t=1:T
    
    %Use the function "update_distribution.m" to construct the t+1
    %distribution using the t one (see the function for details).
    D_new = update_distribution(a_grid,Pr,a_polaa(:,:,t),Dt(:,:,t));
    
    %Construct period t market-clearing condition
    market_cleart(t) = sum(sum(Dt(:,:,t).*a_polaa(:,:,t)));
    %Update the wealth distribution
    Dt(:,:,t+1) = D_new;
    
end

%Return the sequence of market-clearing conditions
y = market_cleart;


