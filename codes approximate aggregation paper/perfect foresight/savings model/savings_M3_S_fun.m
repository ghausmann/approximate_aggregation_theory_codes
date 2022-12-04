%This function takes as inputs the 6 coefficients that characterize the
%M(3) approximation {lV,mV,nV,lM3,mM3,nM3} (responses of changes to the
%aggregate endogenous moments V and E3) to construct a system of 6 independent
%equations.
% The M(3) approximation for the endogenous aggregate states takes the form
%(q-q_ss)   = lV*(V-Vss) + lE3*(3E-E3ss)
%(Vp-Vss)   = mV*(V-Vss) + mM3*(3E-E3ss)
%(M3p-M3ss) = nV*(V-Vss) + nM3*(3E-E3ss)

%where q is the equilibrium bond price, V and M3 are the (pre-determined)
%variance and "skewness" of the wealth distribution, and {q_ss,Vss,M3ss}
%their respective SS values. 
%The inputs are Tmax (maximum number of periods), idgds (average weighted
%derivatives w.r.t. current individual states), idgdq (average derivatives
%w.r.t. future bond prices), and x (vector of coefficients).


function y = savings_M3_S_fun(Tmax,idgds,idgdq,x)

%coefficients q
lV = x(1);
lM3 = x(2);
%coefficients Vp
mV = x(3);
mM3 = x(4);
%coefficients E3p
nV = x(5);
nM3 = x(6);

%Step size for numerical differentiation
h=0.00001;

%Derivatives of tomorrow's moments w.r.t. future asset prices.
idg1dq = idgdq(1,:);
idg2dq = idgdq(2,:);
idg3dq = idgdq(3,:);

%First we calculate the sequence of derivatives of future moments w.r.t
%today's variance and "skewness" (they could be computed by hand, but much
%easier this way..). We then can combine them with the sequences of average
%derivatives w.r.t. future bond prices to compute the total derivative of
%a given tomorrow's moment w.r.t. a given pre-determined one.

%--------------------------------------------------------------------------
%V shock
%--------------------------------------------------------------------------
Vt = zeros(1,Tmax);
M3t = zeros(1,Tmax);
%Small shock to today's variance
Vt(1) = h;

for t = 2:Tmax
    
    Vt(t) = mV*Vt(t-1) + mM3*M3t(t-1);
    M3t(t) =nV*Vt(t-1) + nM3*M3t(t-1);
    
end

%Compute sequence of derivatives of future moments w.r.t today's variance
dVdV_1 = Vt/h;
dM3dV_1 = M3t/h;

%--------------------------------------------------------------------------
%M3 shock
%--------------------------------------------------------------------------
Vt = zeros(1,Tmax);
M3t = zeros(1,Tmax);
%Small shock to today's "skewness"
M3t(1) = h;

for t = 2:Tmax
    
    Vt(t) = mV*Vt(t-1) + mM3*M3t(t-1);
    M3t(t) =nV*Vt(t-1) + nM3*M3t(t-1);
    
end

%Compute sequence of derivatives of future variances and
%covariances w.r.t today's "skewness"
dVdM3_1 = Vt/h;
dM3dM3_1 = M3t/h;
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

%In the paper, I show how to use pertubation to obtain 6 conditions to
%solve for {lV,mV,nV,lE3,mE3,nE3}. The general rule is that the total
%derivative of a given tomorrow's moment of the wealth distribution w.r.t.
%a today's moment is the sum of an average direct effect (as individual
%savings rules depend directly on the moments), and an average indirect
%effect (a small change to a given moment implies a small shock to all the
%individual states, which we must take into account).

%Average indirect effects: of changes to today's Variance and "skewness"
%(it operates through changes to the perturbation shock needed to change
%these moments). See the function "savings_idgds.m" for details.
dBp_dV_ind = idgds(1,1);
dVp_dV_ind = idgds(1,2);
dM3p_dV_ind = idgds(1,3);

dBp_dM3_ind = idgds(2,1);
dVp_dM3_ind = idgds(2,2);
dM3p_dM3_ind = idgds(2,3);


%Average direct effects: total derivatives of tomorrow's moments w.r.t today's
%variance and covariance.
dBp_dV_dir = sum(idg1dq.*(lV*dVdV_1 + lM3*dM3dV_1));
dVp_dV_dir = sum(idg2dq.*(lV*dVdV_1 + lM3*dM3dV_1));
dM3p_dV_dir = sum(idg3dq.*(lV*dVdV_1 + lM3*dM3dV_1));

dBp_dM3_dir = sum(idg1dq.*(lV*dVdM3_1 + lM3*dM3dM3_1));
dVp_dM3_dir = sum(idg2dq.*(lV*dVdM3_1 + lM3*dM3dM3_1));
dE3p_dM3_dir = sum(idg3dq.*(lV*dVdM3_1 + lM3*dM3dM3_1));

%Compute total derivatives
dBpdV = dBp_dV_dir + dBp_dV_ind;
dBpdM3 = dBp_dM3_dir + dBp_dM3_ind;

dVpdV = dVp_dV_dir + dVp_dV_ind;
dVpdM3 = dVp_dM3_dir + dVp_dM3_ind;

dM3pdV = dM3p_dV_dir + dM3p_dV_ind;
dM3pdM3 = dE3p_dM3_dir + dM3p_dM3_ind;

%Since we have 3 aggregate tomorrow's moments (mean, variance and
%"skewness"), and two relevant  pre-determined moments (today's variance
%and "skewness"), we have 6 total derivatives that depend upon
%{lV,mV,nV,lM3,mM3,nM3}, and that in equilibrium must satisfy certain
%conditions.

%Construct the equations
Eq1 = dBpdV;
Eq2 = dBpdM3;
Eq3 = mV - dVpdV;
Eq4 = mM3 - dVpdM3;
Eq5 = nV - dM3pdV;
Eq6 = nM3 - dM3pdM3;

y = [Eq1;Eq2;Eq3;Eq4;Eq5;Eq6];