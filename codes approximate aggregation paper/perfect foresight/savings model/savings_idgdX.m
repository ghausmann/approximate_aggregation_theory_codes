%This function computes a sequence of average derivatives of functions of
%today's savings rule with respect to the sequences of  bond prices and
%borrowing limits, from today until T.
%Inputs are P (vector of model parameters), q (SS bond price), b_grid and
%ui (grids for asset and productivity values with lengths n_b and n_u), Pr
%(the probability transition matrix with size n_u-by-n_u), c_pol (SS
%consumption rule with size n_u-by-n_b), Dss (the stationary distribution
%with size n_u-by-n_b), dS (the step size for numerical differentiation),
%and T (maximum number of periods).

function y = savings_idgdX(P,q,b_grid,ui,Pr,c_pol,Dss,dS,T)

%Model parameters
phi = P(3);

%Sequences of average derivatives wrt. borrowing limits
idgdZ = zeros(1,T); %for the mean
idg2dZ = zeros(1,T); %for the variance
idg3dZ = zeros(1,T); %for the skewness

%Sequences of average derivatives wrt. bond prices
idgdq = zeros(1,T); %for the mean
idg2dq = zeros(1,T); %for the variance
idg3dq = zeros(1,T); %for the skewness

%Sequence of inputs
phi_t = phi*ones(1,T); %borrowing limit
qt = q*ones(1,T); %bond price

%Add the small shocks
phi_t(1) = phi + dS;
qt(1) = q + dS;

c_polp = c_pol; 
c_polq = c_pol; 

b_poli = (exp(ui') + b_grid - c_pol)/q; %SS rule for savings

%Backward iteration to compute derivatives of today's functions w.r.t.
%{phi_0, phi_1, phi_2, ...} and {q0, q1, q2, ...}
for t=1:T
    
    %BORROWING LIMIT
    %Find today's consumption rule given tomorrow's one, and today's
    %borrowing limit.
    c_polip = c_poli_savings_update(P,b_grid,ui,Pr,c_polp,phi_t(t),q);
    %Compute the associated savings rule.
    b_polip = (exp(ui') + b_grid - c_polip)/q;

    %Compute the average derivative of three functions of the savings
    %rule: 
    idgdZ(t) = sum(sum(Dss.*(b_polip-b_poli)/dS));
    idg2dZ(t) = sum(sum(Dss.*(b_polip.^2-b_poli.^2)/dS));
    idg3dZ(t) = sum(sum(Dss.*(b_polip.^3-b_poli.^3)/dS));

    %Update consumption rule
    c_polp=c_polip;
    
    %BOND PRICE
    %Same tasks, but with shocks to bond prices.
    c_poliq = c_poli_savings_update(P,b_grid,ui,Pr,c_polq,phi,qt(t));
    %Compute the associated savings rule.
    b_poliq = (exp(ui') + b_grid - c_poliq)/qt(t);

    %Compute the average derivative of three functions of the savings
    %rule: the savings rule itself (for the mean), (b_poli)^2 (for the
    %variance), and  (b_poli)^3 (for the skewness).
    idgdq(t) = sum(sum(Dss.*(b_poliq-b_poli)/dS));
    idg2dq(t) = sum(sum(Dss.*(b_poliq.^2-b_poli.^2)/dS));
    idg3dq(t) = sum(sum(Dss.*(b_poliq.^3-b_poli.^3)/dS));

    %Update consumption rule
    c_polq=c_poliq;

end

%Return the sequences of average derivatives
y = [idgdZ;idg2dZ;idg3dZ;idgdq;idg2dq;idg3dq];