%This function takes as an input a given value lK from the implied law of
%motion Kp = K_ss + lK1*(K-Kss) + lZ*(Z-1) assumed by households, and
%constructs the corresponding matrix equation as derived in the paper.
%Inputs are P (vector of coefficients), r (SS rental price), T (maximum
%number of periods), idgdX (average derivatives w.r.t. future inputs), and
%lK.

function y = ks_M1_S_fun(P,r,T,idgdX,lK)

%Model parameters
alpha = P(3);

K = (alpha/r)^(1/(1-alpha)); %SS aggregate capital
w = (1-alpha)*(K^alpha)*((1-alpha)); %SS wage rate
%SS partial derivatives of prices w.r.t. aggregate capital
rK = -(1-alpha)*r/K;
wK = alpha*w/K;

%Sequences of average partial derivatives of the savings rule w.r.t. future
%prices
idgdr = idgdX(1,:); %w.r.t the rental price
idgdw = idgdX(2,:); %w.r.t. the wage rate
idgdK = idgdX(3,:); %w.r.t. the wage rate

mytime = 0:1:(T-1);

%Explicit equation (see the paper)
dYp_dK = sum((rK*idgdr + wK*idgdw + idgdK).*(lK.^mytime));

y = dYp_dK;