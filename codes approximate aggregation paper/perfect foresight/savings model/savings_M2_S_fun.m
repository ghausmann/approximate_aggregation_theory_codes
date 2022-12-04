%This function takes as inputs the 2 coefficients that characterize the
%M(2) approximation {lV,mV} (responses of changes to the variance) to
%construct a system of 2 independent equations.
% The M(2) approximation for the variance takes the form
%(q-q_ss) = lV*(V-Vss)
%(Vp-Vss) = mV*(V-Vss)
%where q is the equilibrium bond price, V the (pre-determined) variance of
%the wealth distribution, and {q_ss,Vss} their respective SS values.
%Inputs are Tmax (maximum number of periods), idGds (average weighted
%derivatives w.r.t. current individual states), idGdq (average derivatives
%w.r.t. future bond prices), and x (vector of coefficients {lV,mV}).

function y = savings_M2_S_fun(Tmax,idGds,idgdq,x)


lV = x(1); %coefficient q
mV = x(2); %coefficient V

%Derivatives of tomorrow's moments w.r.t. future bond prices.
idg1dq = idgdq(1,:); %For the mean
idg2dq = idgdq(2,:); %For the variance

mytime = 0:1:(Tmax-1);

%--------------------------------------------------------------------------
%Explicit formulas for the total direct effects (see the paper)
%--------------------------------------------------------------------------
dEb_dV_dir = lV*sum(idg1dq.*(mV.^mytime));
dVt_dV_dir = lV*sum(idg2dq.*(mV.^mytime));


dEb_dV_ind = idGds(1);
dVt_dV_ind = idGds(2);

%Total derivatives of tomorrow's mean and variance w.r.t. today's variance.
dEb_dV_total = dEb_dV_dir + dEb_dV_ind;
dVt_dV_total = dVt_dV_dir + dVt_dV_ind;

%Form the equations
Eq1 = dEb_dV_total;
Eq2 = mV - dVt_dV_total;

y = [Eq1;Eq2];