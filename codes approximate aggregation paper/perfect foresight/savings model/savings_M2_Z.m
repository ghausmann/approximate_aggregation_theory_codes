%This function takes as inputs the 2 coefficients that characterize the
%M(2) approximation {lZ,mZ} (responses to exogenous states) to construct a
%system of 2 independent equations.
% The full M(2) approximation takes the form
%(q-q_ss) = lV*(V-Vss) + lZ*(Z-Zss)
%(Vp-Vss) = mV*(V-Vss) + mZ*(Z-Zss)
%where q is the equilibrium bond price, V the (pre-determined) variance of
%the wealth distribution, {q_ss,Vss} their respective SS values, and Z the
%exogenous state. Here we take as given the coefficients {lV,mV} (responses
%of changes to the aggregate endogenous moments), as they can be solved
%independently and in advance.
%Inputs are rho (persistence of the exogenous process), v (vector of coefficients
%{lV,mV}), Tmax (maximum number of periods), idGdq (average
%derivatives w.r.t. future bond prices), and idGdZ (average derivatives w.r.t.
%future exogenous shocks).

function y = savings_M2_Z(rho,v,Tmax,idgdq,idgdZ)

%Already solved coefficients
lV = v(1);
mV = v(2);

%Step size for numerical differentiation
d=0.00001;

%Derivatives of tomorrow's moments w.r.t. future bond prices.
idg1dq = idgdq(1,:);
idg2dq = idgdq(2,:);

%Derivatives of tomorrow's moments w.r.t. future exogenous shocks.
idg1dZ = idgdZ(1,:);
idg2dZ = idgdZ(2,:);

%First we calculate the sequence of derivatives of future moments and
%exogenous state values w.r.t. today's exogenous shock (they could be
%computed by hand, but much easier this way..).
Vt = zeros(1,Tmax);
Zt = zeros(1,Tmax);
Zt(1) = d;

%Fix values for coefficient equal to one
mZ = 1;

for t = 2:Tmax
    
    Zt(t) =  rho*Zt(t-1);    
    Vt(t) = mV*Vt(t-1) +  mZ*Zt(t-1);
        
end

dVdZ = Vt/d;
dbldZ = Zt/d;

Bp0 = sum(idg1dZ.*dbldZ);
Vp0 = sum(idg2dZ.*dbldZ);

Bpb = sum(idg1dq.*(dbldZ));
Vpb = sum(idg2dq.*(dbldZ));

%These objects are such that they are proportional to mZ. Since we fixed
%mZ=1, we already have the desired values.
Bpv = lV*sum(idg1dq.*dVdZ);
Vpv = lV*sum(idg2dq.*dVdZ);

%We have all the objects needed to solve for the coefficients with
%straightforward matrix algebra. 
Mc = [Bpb Bpv;Vpb (Vpv-1)];
Mf = -[Bp0;Vp0];

msol = (Mc^-1)*Mf;

%Return the coefficients
y = msol';