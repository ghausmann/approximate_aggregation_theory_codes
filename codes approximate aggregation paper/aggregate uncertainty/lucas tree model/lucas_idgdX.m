%This functions computes a sequence of average derivatives of today's asset
%rule with respect to sequences of output and the asset price from today
%until T).
%Inputs are P (vector of model parameters), pss (SS asset price), a_grid
%and zi (grids for asset and productivity values with lengths n_a and n_z),
%Pr (the probability transition matrix with size n_z-by-n_z), c_pol (SS
%consumption rule with size n_z-by-n_a), a_poli (SS savings rule with size
%n_z-by-n_a), Dss (the stationary distribution with size n_z-by-n_a), dS
%(the step size for numerical differentiation), and T (maximum number of
%periods).

function y = lucas_idgdX(P,pss,a_grid,zi,Pr,c_pol,Dss,dS,T)

%Model parameters
alpha = P(3);

%Sequences of average derivatives
idgdY = zeros(1,T);
idgdp = zeros(1,T);

%Sequences of output and prices
Yt = ones(1,T);
Yt1 = ones(1,T); %dYt1 is used for a shock to tomorrow's output
pt = pss*ones(1,T);
pt1 = pss*ones(1,T); %pt1 is used for a shock to tomorrow's price

Yt(1) = 1 + dS;
Yt1(2) = 1 +dS;
pt(1) = pss + dS;
pt1(2) = pss + dS;

c_polY = c_pol;
c_polp = c_pol;

%SS asset share rule
a_poli = ( (1+zi')*(1-alpha) + (1+a_grid)*(pss+alpha) - c_pol)/pss - 1;

%Backward iteration to compute derivatives of today's savings rule w.r.t.
%{Y0, Y1, Y2, ....} and {p0, p1, p2, ....}.
for n = 1:T
    
    %OUTPUT
    %Find today's consumption rule given tomorrow's one, and today's and
    %tomorrow's output.
    c_poliY =c_poli_lucas_update(P,a_grid,zi,Pr,c_polY,Yt(n),Yt1(n),pss,pss);
    %Compute the associated asset rule.
    a_polid = ( (1+zi')*(1-alpha)*(Yt(n) ) + (1+a_grid)*((pss )+alpha*(Yt(n))) - c_poliY)/(pss ) - 1;
    %Compute the average derivative of of the asset rule
    idgdY(n) = sum(sum(Dss.*(a_polid-a_poli)/dS));
    %Update consumption rule
    c_polY = c_poliY;
    
    %PRICE
    %Same tasks, but with shocks to prices
    c_polip =c_poli_lucas_update(P,a_grid,zi,Pr,c_polp,1,1,pt(n),pt1(n));
    a_polip = ( (1+zi')*(1-alpha) + (1+a_grid)*((pt(n) )+alpha) - c_polip)/(pt(n) ) - 1;
    idgdp(n) = sum(sum(Dss.*(a_polip-a_poli)/dS));
    c_polp = c_polip;
    
end

%Return the sequences of average derivatives
y = [idgdY;idgdp];


