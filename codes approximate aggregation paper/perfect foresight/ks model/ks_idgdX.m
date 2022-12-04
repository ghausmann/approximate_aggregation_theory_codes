%This function computes a sequence of average derivatives of today's asset
%share rule with respect to sequences of aggregate inputs (prices and
%aggregate capital) from today until terminal period T.
%Inputs are P (vector of model parameters), r (SS rental price), a_grid and
%zi (grids for asset and productivity values with lengths n_a and n_z), Pr
%(probability transition matrix with size nz-by-nz), c_pol (ss consumption
%rule with size n_z-by-n_a), Dss (the stationary distribution  with size
%n_z-by-n_a), dS (the step size for numerical differentiation), and T
%(maximum number of periods).

function y = ks_idgdX(P,r,a_grid,zi,Pr,c_pol,Dss,dS,T)

%Model parameters
alpha = P(3);
d = P(4);

K = (alpha/r)^(1/(1-alpha)); %SS aggregate capital
w = (1-alpha)*(K^alpha)*((1-alpha)); %SS wage rate

%Sequences of average derivatives
idgdr = zeros(1,T); %w.r.t. future rental prices
idgdw = zeros(1,T); %w.r.t. future wage rates
idgdK = zeros(1,T); %w.r.t. future aggregate capital

%Sequences of prices and capital
rt = r*ones(1,T);
rt1 = r*ones(1,T); %rt1 is used for a shock to tomorrow's price
wt = w*ones(1,T);
Kt = K*ones(1,T);
Kt1 = K*ones(1,T); %Kt1 is used for a shock to tomorrow's capital

%Add the small shocks
rt(1) = r + dS;
wt(1) = w + dS;
rt1(2) = r + dS;
Kt(1) = K + dS;
Kt1(2) = K + dS;

c_polr = c_pol;
c_polw = c_pol;
c_polK = c_pol;

%SS asset share rule
a_poli = ((1+zi)'*w + K*(1+r-d)*(a_grid+1) - c_pol)/K - 1;

%Backward iteration to compute derivatives of today's savings rule w.r.t.
%{r0, r1, r2, ...}, {w0,w1, w2, ...}, and {K0, K1, K2, ...}.
for n = 1:T

    %RENTAL PRICES
    %Find today's consumption rule given tomorrow's one, and today's and
    %tomorrow's prices and capital.
    c_polir =c_poli_ks_update(P,a_grid,zi,Pr,c_polr,rt(n),rt1(n),w,K,K);
    %Compute the associated asset rule.
    a_polir = ((w  )*(1+zi') + K*(1+(rt(n))-d)*(a_grid+1) - c_polir)/K - 1;
    %Compute the average derivative of the savings rule
    idgdr(n) = sum(sum(Dss.*(a_polir-a_poli)/dS));
    %Update consumption rule
    c_polr = c_polir;

    %WAGE RATES
    %Same tasks, but with shocks to wages.
    c_poliw =c_poli_ks_update(P,a_grid,zi,Pr,c_polw,r,r,wt(n),K,K);
    a_poliw = ((wt(n) )*(1+zi') + K*(1+(r )-d)*(a_grid+1) - c_poliw)/K - 1;
    idgdw(n) = sum(sum(Dss.*(a_poliw-a_poli)/dS));
    c_polw = c_poliw;

    %AGGREGATE CAPITAL
    %Same tasks, but with shocks to aggregate capital stock.
    c_poliK =c_poli_ks_update(P,a_grid,zi,Pr,c_polK,r,r,w,Kt(n),Kt1(n));
    a_poliK = ((w  )*(1+zi') + Kt(n)*(1+r-d)*(a_grid+1) - c_poliK)/Kt1(n) - 1;
    idgdK(n) = sum(sum(Dss.*(a_poliK-a_poli)/dS));
    c_polK = c_poliK;

end

%Return the sequences of average derivatives
y = [idgdr;idgdw;idgdK];