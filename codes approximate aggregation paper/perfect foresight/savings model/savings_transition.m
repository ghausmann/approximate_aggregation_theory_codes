%This function computes a sequence of size T of excess demands for each
%period bond market, given sequences of candidate bond prices and borrowing
%limits. The equilibrium sequence of bond prices is the one that clears all
%markets (makes all the T excess demands 0).
%Inputs are P (vector of model parameters), T (number of periods), phi_t
%(sequence of borrowing limits), the grids b_grid (with length n_b) and ui
%(with length S), Pr (Probability transition matrix with size n_u-by-n_u),
%c_pol (the SS policy rule for consumption of size n_u-by-n_b), the
%stationary wealth distribution D0 with size n_u-by-n_b, the SS bond price,
%and the sequence of bond prices x.

function y  = savings_transition(P,T,phi_t,b_grid,ui,Pr,c_pol,D0,qss_ha,x)

%Model parameters
phi_ss = P(3);

n_u = length(ui);
n_b = length(b_grid);

%Sequence of bond prices
qt = x;

%--------------------------------------------------------------------------
%BACKWARD ITERATION (of policy rules)
%--------------------------------------------------------------------------

%Compute the sequence of policy functions from the Terminal period T to the
%current period t0, using backward iteration

%SS saving rule
b_poli = ( exp(ui') + b_grid - c_pol)/qss_ha;
b_poli(b_poli<=-phi_ss) = -phi_ss;
b_poli(b_poli>=b_grid(n_b)) = b_grid(n_b);

b_polaa = zeros(n_u,n_b,T+1);
b_polaa(:,:,T+1) = b_poli;

c_polD = c_pol;

for n=(T):-1:1

    %Use the function "c_poli_savings_update.m" to iterate backwards (based
    %on the endogenous grid method)

    c_poli_dD = c_poli_savings_update(P,b_grid,ui,Pr,c_polD,phi_t(n),qt(n));
    %Compute the correspoding savings rule
    b_poli_dD = (exp(ui') + b_grid - c_poli_dD)/qt(n);
    b_poli_dD(b_poli_dD<=(-phi_t(n))) = -phi_t(n);
    b_poli_dD(b_poli_dD>=b_grid(n_b)) = b_grid(n_b);

    c_polD = c_poli_dD;
    b_polaa(:,:,n) = b_poli_dD;
end

%--------------------------------------------------------------------------
%FORWARD ITERATION (of wealth distribution)
%--------------------------------------------------------------------------

%Sequence of wealth distributions
Dt = zeros(n_u,n_b,(T+1));
%Sequence of market-clearing conditions
market_cleart = zeros(1,(T));
Dt(:,:,1) = D0;

%Iterate forward
for t=1:T

    %Use the function "update_distribution.m" to construct the t+1
    %distribution using the t one (see the function for details).
    D_new = update_distribution(b_grid,Pr,b_polaa(:,:,t),Dt(:,:,t));

    %Construct period t market-clearing condition
    market_cleart(t) = sum(sum(Dt(:,:,t).*b_polaa(:,:,t)));
    %Update the wealth distribution
    Dt(:,:,t+1) = D_new;

end

%Return the sequence of market-clearing conditions
y = market_cleart;


