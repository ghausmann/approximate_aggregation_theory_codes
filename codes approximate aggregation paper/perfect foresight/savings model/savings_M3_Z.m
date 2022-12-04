%This function solves for the 3 coefficients that characterize the
%M(3) approximation {lZ,mZ,nZ} using basic matrix algebra.
% The full M(3) approximation takes the form
%(q-q_ss)   = lV*(V-Vss) + lM3*(E3-E3ss) + lZ*(Z-Zss)
%(Vp-Vss)   = mV*(V-Vss) + mM3*(E3-E3ss) + mZ*(Z-Zss)
%(M3p-M3ss) = nV*(V-Vss) + nM3*(E3-E3ss) + nZ*(Z-Zss)
%where q is the equilibrium bond price, V and M3 are the (pre-determined)
%variance and "skewness" of the wealth distribution, Z the exogenous state,
%and {q_ss,Vss,M3ss,Zss} their respective SS values. Here we take as given
%the coefficients {lV,mV,nV,lM3,mM3,nM3} (responses of changes to the
%aggregate endogenous moments), as they can be solved independently and in
%advance.
%Inputs are rho (persistence of the exogenous process), v (the vector of
%coefficients), Tmax (maximum number of periods), idgdq (average
%derivatives w.r.t. future bond prices), and idgdZ (average derivatives
%w.r.t. future exogenous shocks).

function y = savings_M3_Z(rho,v,Tmax,idgdq,idgdZ)

%Already solved coefficients
lV = v(1);
lM3 = v(2);
mV = v(3);
mM3 = v(4);
nV = v(5);
nM3 = v(6);

%Step size for numerical differentiation
d=0.00001;

%Derivatives of tomorrow's moments w.r.t. future bond prices.
idg1dq = idgdq(1,:);
idg2dq = idgdq(2,:);
idg3dq = idgdq(3,:);

%Derivatives of tomorrow's moments w.r.t. future exogenous shocks.
idg1dZ = idgdZ(1,:);
idg2dZ = idgdZ(2,:);
idg3dZ = idgdZ(3,:);

%First we calculate the sequence of derivatives of future moments and
%exogenous state values w.r.t. today's exogenous shock (they could be
%computed by hand, but much easier this way..).
Vt = zeros(1,Tmax);
M3t = zeros(1,Tmax);
Zt = zeros(1,Tmax);
Zt(1) = d;

%Fix values for coefficients
mZ = 1;
nZ = 2;
for t = 2:Tmax

    Zt(t) = rho*Zt(t-1);
    Vt(t) = mV*Vt(t-1) + mM3*M3t(t-1) +  mZ*Zt(t-1);
    M3t(t) =nV*Vt(t-1) + nM3*M3t(t-1) + nZ*Zt(t-1);

end

dVdZ0 = Vt/d;
dM3dZ0 = M3t/d;
dZdZ0 = Zt/d;

%These objects are such that they can be decomposed in two components, one
%proportional to mZ, and the other to nZ.
obj1 = sum(idg1dq.*(lV*dVdZ0 + lM3*dM3dZ0 ));
obj2 = sum(idg2dq.*(lV*dVdZ0 + lM3*dM3dZ0));
obj3 = sum(idg3dq.*(lV*dVdZ0 + lM3*dM3dZ0));

%Repeat the same process, but with different values for mZ and nZ.
Vt = zeros(1,Tmax);
M3t = zeros(1,Tmax);

mZ = 3;
nZ = 4;
for t = 2:Tmax


    Vt(t) =mV*Vt(t-1) + mM3*M3t(t-1) + mZ*Zt(t-1);
    M3t(t)=nV*Vt(t-1) + nM3*M3t(t-1) + nZ*Zt(t-1);

end

dVdZ0 = Vt/d;
dM3dZ0 = M3t/d;

%Compute the new objects that can be decomposed in two components, the same
%ones as before!
obj11 = sum(idg1dq.*(lV*dVdZ0 + lM3*dM3dZ0 ));
obj22 = sum(idg2dq.*(lV*dVdZ0 + lM3*dM3dZ0));
obj33 = sum(idg3dq.*(lV*dVdZ0 + lM3*dM3dZ0));

%Now we have identification, so solve for the components proportional to mZ
%and nZ
comps1 = ([1 2;3 4]^-1)*[obj1;obj11];
comps2 = ([1 2;3 4]^-1)*[obj2;obj22];
comps3 = ([1 2;3 4]^-1)*[obj3;obj33];

%We have all the objects needed to solve for the coefficients with
%straightforward matrix algebra.
Bp0 = sum(idg1dZ.*dZdZ0);
Vp0 = sum(idg2dZ.*dZdZ0);
Mp0 = sum(idg3dZ.*dZdZ0);

Bpb = sum(idg1dq.*(dZdZ0));
Vpb = sum(idg2dq.*(dZdZ0));
Mpb = sum(idg3dq.*(dZdZ0));

Bpv = comps1(1);
Vpv = comps2(1);
Mpv = comps3(1);

Bpe = comps1(2);
Vpe = comps2(2);
Mpe = comps3(2);


Cc = [Bpb Bpv Bpe;Vpb (Vpv-1) Vpe;Mpb Mpv (Mpe-1)];
Cf = -[Bp0;Vp0;Mp0];

msol = (Cc^-1)*Cf;

%Return the coefficients
y = msol';