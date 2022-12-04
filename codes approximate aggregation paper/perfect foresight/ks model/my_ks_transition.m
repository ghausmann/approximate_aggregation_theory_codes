%This function computes a sequence of size T of excess demands (one for
%each period), given a sequence of candidate rental prices and the initial
%aggregate productivity shock. The equilibrium sequence of prices is the
%one that clears all markets.
%The inputs are P (vector of model parameters), T (number of periods), dD
%(initial productivity shock), the grids a_grid and zi (with lengths n_a
%and n_z), Pr (Probability transition matrix with size n_z-by-n_z), c_pol
%(the SS policy rule for consumption with size n_z-by-n_a), rss (SS rental
%price), Dss (the stationary wealth distribution with size n_z-by-n_a),
%and x, a vector of candidate rental prices.

function y  = my_ks_transition(P,T,dD,a_grid,zi,Pr,c_pol,rss,Dss,x)

%Model parameters
alpha = P(3);
d = P(4);
rho = P(5);

%SS capital and wage rate
Kss = (alpha/rss)^(1/(1-alpha));
wss = ((1-alpha))*(1-alpha)*Kss^alpha;

%Today's rental price
r0 = alpha*(1+dD)*(Kss)^(alpha-1);

n_a = length(zi);
n_z = length(a_grid);

%Sequence of candidate rental prices
%Note: r0 and rss already known!
rt = [r0 x rss];
%Excess productivity
dZ = zeros(1,T+1);
dZ(1) = dD;

%Compute sequence for aggregate productivity
for t=2:T+1

    dZ(t) = rho*dZ(t-1);

end
Zt = 1 + dZ;

%Compute associated sequences for aggregate capital and the wage rate
Kt_1 = (alpha*Zt./rt).^(1/(1-alpha));
wt = (1-alpha)*Zt.*(Kt_1.^alpha)*((1-alpha));

%--------------------------------------------------------------------------
%BACKWARD ITERATION (of policy rules)
%--------------------------------------------------------------------------

%SS policy rule for asset shares
a_poli = ((1+zi')*wss + Kss*(1+rss-d)*(a_grid+1) - c_pol)/Kss - 1;
%Make sure the constraints hold.
a_poli(a_poli<=-1) = -1;
a_poli(a_poli>=(a_grid(n_z))) = a_grid(n_z);

%Compute the sequence of policy functions from the Terminal period T+1 to the
%current period t0, using backward iteration
c_pol_dD  = c_pol;

a_polaa = zeros(n_a,n_z,T+1);
a_polaa(:,:,T+1) = a_poli;

for n=(T):-1:1
    %Use the function "c_poli_ks_update.m" to iterate backwards.
    c_poli_dD =  c_poli_ks_update(P,a_grid,zi,Pr,c_pol_dD,rt(n),rt(n+1),wt(n),Kt_1(n),Kt_1(n+1));
    %Compute the correspoding asset rule
    a_poli_dD = (wt(n)*(1+zi') + Kt_1(n)*(1+rt(n)-d)*(a_grid+1) - c_poli_dD)/Kt_1(n+1) - 1;
    a_poli_dD(a_poli_dD<=-1) = -1;
    a_poli_dD(a_poli_dD>=(a_grid(n_z))) = a_grid(n_z);

    c_pol_dD= c_poli_dD; %update the consumption rule
    a_polaa(:,:,n) = a_poli_dD; %store the asset rule
end

%--------------------------------------------------------------------------
%FORWARD ITERATION (of wealth distribution)
%--------------------------------------------------------------------------

market_cleart = zeros(1,T-1); %Sequence of market-clearing conditions
D0 = Dss; %Initial wealth distribution (the SS one).

%Iterate forward
for t=1:T-1

    %Construct period t market-clearing condition (excess demand)
    Mp = sum(sum(D0.*a_polaa(:,:,t)));
    market_cleart(t) = Mp;

    %Update the wealth distribution (see the function
    %"update_distribution.m" for details).
    D1 = update_distribution(a_grid,Pr,a_polaa(:,:,t),D0);
    D0 = D1;
end

%Return the sequence of market-clearing conditions
y = market_cleart;