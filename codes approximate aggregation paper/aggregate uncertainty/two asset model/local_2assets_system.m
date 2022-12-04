%This function evaluates the system of equations that solves the two-asset
%model.
%The inputs are  P (vector of parameters), x (vector of coefficients), and
%ut (vector of T simulated idiosyncratic innovations).
%The output are the equations evaluated at the inputs.

function y = local_2assets_system(P,x,ut)

%parameters
beta = P(1);
gamma = P(2);
phi = P(3);
su = P(4);
sz = P(5);
rz = P(6);

%coefficients
qss = x(1);
nz = x(2);
pss = x(3);
lz = x(4);
gb = x(5);
ga = x(6);
gu = x(7);
hb = x(8);
ha = x(9);
hu = x(10);
gss = x(11);
hss = x(12);

%coefficients of the discriminant function
fss = pss*gss + qss*hss;
fb = pss*gb + qss*hb;
fa = pss*ga + qss*ha;
fu = pss*gu + qss*hu;


%% 
%--------------------------------------------------------------------------
%Simulation part
%--------------------------------------------------------------------------

T = length(ut);
bt = zeros(1,T); %safe asset
at = zeros(1,T); %risky asset

at(1) = gss;
bt(1) = hss;

%Simulate a large time-series to approximate the stationary distribution
for t=2:T
    
    at(t) = gss + gb*(bt(t-1) - hss) + ga*(at(t-1) - gss) + gu*ut(t);
    bt(t) = hss + hb*(bt(t-1) - hss) + ha*(at(t-1) - gss) + hu*ut(t);
    
    %Correct values if the constrain binds
    ft = qss*bt(t) + pss*at(t);
    if ft<=-phi
        ub = -( phi + fss + fb*(bt(t-1)-hss) + fa*(at(t-1)-gss) )/fu;
        
        at(t) = gss + gb*(bt(t-1) - hss) + ga*(at(t-1) - gss)  + gu*ub;
        bt(t) = hss + hb*(bt(t-1) - hss) + ha*(at(t-1) - gss)  + hu*ub;
        
    end
    
end

at = at(1001:T);
bt = bt(1001:T);

%Market-clearing equations
ZOb = mean(bt);
ZOa = mean(at);

%% 
%--------------------------------------------------------------------------
%Probability objects
%--------------------------------------------------------------------------

%Compute the moments of the truncated normal distribution. See Greene
%(Econometric analysis), Chapter 19th for details.

a = -(phi+fss)/(fu*su );

phib = 0.5*(1 + erf((a)/(2^0.5)));
phil = (1/((2*pi)^0.5))*exp(-0.5*(( (a))^2));
lau = phil/(1-phib);
lac = -phil/phib;
dau = lau*(lau-a);
dac = lac*(lac-a);

Prc = phib; %Probability of constrained tomorrow

u1u = lau; %mean of the individual shock if unconstrained
u1c = lac; %mean of the individual shock if constrained

vu_u = 1 - dau; %variance of the individual shock if unconstrained
vu_c = 1 - dac; %variance of the individual shock if constrained

%% 

%--------------------------------------------------------------------------
%Conditions from the Euler equations
%--------------------------------------------------------------------------

U1 =(Prc*beta)/(gss + hss + phi + su*u1c + 1)^gamma - qss/(gss + hss - gss*pss - hss*qss + 1)^gamma - (beta*(Prc - 1))/(gss + hss - pss*(gss + gu*su*u1u) + su*u1u - qss*(hss + hu*su*u1u) + 1)^gamma - (beta*gamma*(Prc - 1)*(gamma + 1)*(gss*sz - sz + gss*lz*sz + hss*nz*sz)^2)/(2*(gss + hss - pss*(gss + gu*su*u1u) + su*u1u - qss*(hss + hu*su*u1u) + 1)^(gamma + 2)) + (Prc*beta*gamma*(gamma + 1)*(gss*sz - sz + gss*lz*sz + hss*nz*sz)^2)/(2*(gss + hss + phi + su*u1c + 1)^(gamma + 2)) - (beta*gamma*vu_u*(Prc - 1)*(gamma + 1)*(gu*pss*su - su + hu*qss*su)^2)/(2*(gss + hss - pss*(gss + gu*su*u1u) + su*u1u - qss*(hss + hu*su*u1u) + 1)^(gamma + 2)) + (Prc*beta*gamma*su^2*vu_c*(gamma + 1))/(2*(gss + hss + phi + su*u1c + 1)^(gamma + 2));
U2 =(beta*((Prc*gamma*(gamma + 1)*(gss*sz - sz + gss*lz*sz + hss*nz*sz)^2)/(gss + hss + phi + su*u1c + 1)^(gamma + 2) - (gamma*(Prc - 1)*(gamma + 1)*(gss*sz - sz + gss*lz*sz + hss*nz*sz)^2)/(gss + hss - pss*(gss + gu*su*u1u) + su*u1u - qss*(hss + hu*su*u1u) + 1)^(gamma + 2)))/2 - beta*((Prc - 1)/(gss + hss - pss*(gss + gu*su*u1u) + su*u1u - qss*(hss + hu*su*u1u) + 1)^gamma - Prc/(gss + hss + phi + su*u1c + 1)^gamma) - pss/(gss + hss - gss*pss - hss*qss + 1)^gamma - beta*sz*((Prc*gamma*(gss*sz - sz + gss*lz*sz + hss*nz*sz))/(gss + hss + phi + su*u1c + 1)^(gamma + 1) - (gamma*(Prc - 1)*(gss*sz - sz + gss*lz*sz + hss*nz*sz))/(gss + hss - pss*(gss + gu*su*u1u) + su*u1u - qss*(hss + hu*su*u1u) + 1)^(gamma + 1)) - (beta*gamma*vu_u*(Prc - 1)*(gamma + 1)*(gu*pss*su - su + hu*qss*su)^2)/(2*(gss + hss - pss*(gss + gu*su*u1u) + su*u1u - qss*(hss + hu*su*u1u) + 1)^(gamma + 2)) + (Prc*beta*gamma*su^2*vu_c*(gamma + 1))/(2*(gss + hss + phi + su*u1c + 1)^(gamma + 2));
U1b =(beta*gamma*(Prc - 1)*(gb + hb - qss*(hb^2 + gb*ha) - pss*(ga*gb + gb*hb)))/(gss + hss - pss*(gss + gu*su*u1u) + su*u1u - qss*(hss + hu*su*u1u) + 1)^(gamma + 1) - (Prc*beta*gamma*(gb + hb))/(gss + hss + phi + su*u1c + 1)^(gamma + 1) - (gamma*qss*(gb*pss + hb*qss - 1))/(gss + hss - gss*pss - hss*qss + 1)^(gamma + 1) + (Prc*beta*gamma*gb*sz*(gamma + 1)*(gss*sz - sz + gss*lz*sz + hss*nz*sz))/(gss + hss + phi + su*u1c + 1)^(gamma + 2) - (Prc*beta*gamma*(gb + hb)*(gamma + 1)*(gamma + 2)*(gss*sz - sz + gss*lz*sz + hss*nz*sz)^2)/(2*(gss + hss + phi + su*u1c + 1)^(gamma + 3)) - (beta*gamma*gb*sz*(Prc - 1)*(gamma + 1)*(gss*sz - sz + gss*lz*sz + hss*nz*sz))/(gss + hss - pss*(gss + gu*su*u1u) + su*u1u - qss*(hss + hu*su*u1u) + 1)^(gamma + 2) + (beta*gamma*(Prc - 1)*(gamma + 1)*(gamma + 2)*(gss*sz - sz + gss*lz*sz + hss*nz*sz)^2*(gb + hb - qss*(hb^2 + gb*ha) - pss*(ga*gb + gb*hb)))/(2*(gss + hss - pss*(gss + gu*su*u1u) + su*u1u - qss*(hss + hu*su*u1u) + 1)^(gamma + 3)) + (beta*gamma*vu_u*(Prc - 1)*(gamma + 1)*(gamma + 2)*(gu*pss*su - su + hu*qss*su)^2*(gb + hb - qss*(hb^2 + gb*ha) - pss*(ga*gb + gb*hb)))/(2*(gss + hss - pss*(gss + gu*su*u1u) + su*u1u - qss*(hss + hu*su*u1u) + 1)^(gamma + 3)) - (Prc*beta*gamma*su^2*vu_c*(gb + hb)*(gamma + 1)*(gamma + 2))/(2*(gss + hss + phi + su*u1c + 1)^(gamma + 3));
U1a =(beta*gamma*(Prc - 1)*(ga + ha - pss*(ga^2 + gb*ha) - qss*(ga*ha + ha*hb)))/(gss + hss - pss*(gss + gu*su*u1u) + su*u1u - qss*(hss + hu*su*u1u) + 1)^(gamma + 1) - (Prc*beta*gamma*(ga + ha))/(gss + hss + phi + su*u1c + 1)^(gamma + 1) - (gamma*qss*(ga*pss + ha*qss - 1))/(gss + hss - gss*pss - hss*qss + 1)^(gamma + 1) + (Prc*beta*ga*gamma*sz*(gamma + 1)*(gss*sz - sz + gss*lz*sz + hss*nz*sz))/(gss + hss + phi + su*u1c + 1)^(gamma + 2) - (Prc*beta*gamma*(ga + ha)*(gamma + 1)*(gamma + 2)*(gss*sz - sz + gss*lz*sz + hss*nz*sz)^2)/(2*(gss + hss + phi + su*u1c + 1)^(gamma + 3)) - (beta*ga*gamma*sz*(Prc - 1)*(gamma + 1)*(gss*sz - sz + gss*lz*sz + hss*nz*sz))/(gss + hss - pss*(gss + gu*su*u1u) + su*u1u - qss*(hss + hu*su*u1u) + 1)^(gamma + 2) + (beta*gamma*(Prc - 1)*(gamma + 1)*(gamma + 2)*(gss*sz - sz + gss*lz*sz + hss*nz*sz)^2*(ga + ha - pss*(ga^2 + gb*ha) - qss*(ga*ha + ha*hb)))/(2*(gss + hss - pss*(gss + gu*su*u1u) + su*u1u - qss*(hss + hu*su*u1u) + 1)^(gamma + 3)) + (beta*gamma*vu_u*(Prc - 1)*(gamma + 1)*(gamma + 2)*(gu*pss*su - su + hu*qss*su)^2*(ga + ha - pss*(ga^2 + gb*ha) - qss*(ga*ha + ha*hb)))/(2*(gss + hss - pss*(gss + gu*su*u1u) + su*u1u - qss*(hss + hu*su*u1u) + 1)^(gamma + 3)) - (Prc*beta*gamma*su^2*vu_c*(ga + ha)*(gamma + 1)*(gamma + 2))/(2*(gss + hss + phi + su*u1c + 1)^(gamma + 3));
U1ut =(beta*gamma*(Prc - 1)*(gu + hu - pss*(ga*gu + gb*hu) - qss*(gu*ha + hb*hu)))/(gss + hss - pss*(gss + gu*su*u1u) + su*u1u - qss*(hss + hu*su*u1u) + 1)^(gamma + 1) - (Prc*beta*gamma*(gu + hu))/(gss + hss + phi + su*u1c + 1)^(gamma + 1) - (gamma*qss*(gu*pss + hu*qss - 1))/(gss + hss - gss*pss - hss*qss + 1)^(gamma + 1) + (Prc*beta*gamma*gu*sz*(gamma + 1)*(gss*sz - sz + gss*lz*sz + hss*nz*sz))/(gss + hss + phi + su*u1c + 1)^(gamma + 2) + (beta*gamma*(Prc - 1)*(gamma + 1)*(gamma + 2)*(gu + hu - pss*(ga*gu + gb*hu) - qss*(gu*ha + hb*hu))*(gss*sz - sz + gss*lz*sz + hss*nz*sz)^2)/(2*(gss + hss - pss*(gss + gu*su*u1u) + su*u1u - qss*(hss + hu*su*u1u) + 1)^(gamma + 3)) - (Prc*beta*gamma*(gu + hu)*(gamma + 1)*(gamma + 2)*(gss*sz - sz + gss*lz*sz + hss*nz*sz)^2)/(2*(gss + hss + phi + su*u1c + 1)^(gamma + 3)) - (beta*gamma*gu*sz*(Prc - 1)*(gamma + 1)*(gss*sz - sz + gss*lz*sz + hss*nz*sz))/(gss + hss - pss*(gss + gu*su*u1u) + su*u1u - qss*(hss + hu*su*u1u) + 1)^(gamma + 2) + (beta*gamma*vu_u*(Prc - 1)*(gamma + 1)*(gamma + 2)*(gu*pss*su - su + hu*qss*su)^2*(gu + hu - pss*(ga*gu + gb*hu) - qss*(gu*ha + hb*hu)))/(2*(gss + hss - pss*(gss + gu*su*u1u) + su*u1u - qss*(hss + hu*su*u1u) + 1)^(gamma + 3)) - (Prc*beta*gamma*su^2*vu_c*(gu + hu)*(gamma + 1)*(gamma + 2))/(2*(gss + hss + phi + su*u1c + 1)^(gamma + 3));
U1zt =(Prc*beta*gamma*(gss*rz - rz + gss*lz*rz + hss*nz*rz))/(gss + hss + phi + su*u1c + 1)^(gamma + 1) - (gamma*qss*(gss + gss*lz + hss*nz - 1))/(gss + hss - gss*pss - hss*qss + 1)^(gamma + 1) - (beta*gamma*(Prc - 1)*(gss*rz - rz + gss*lz*rz + hss*nz*rz))/(gss + hss - pss*(gss + gu*su*u1u) + su*u1u - qss*(hss + hu*su*u1u) + 1)^(gamma + 1) - nz/(gss + hss - gss*pss - hss*qss + 1)^gamma + (Prc*beta*gamma*(gamma + 1)*(gamma + 2)*(gss*rz - rz + gss*lz*rz + hss*nz*rz)*(gss*sz - sz + gss*lz*sz + hss*nz*sz)^2)/(2*(gss + hss + phi + su*u1c + 1)^(gamma + 3)) - (beta*gamma*(Prc - 1)*(gamma + 1)*(gamma + 2)*(gss*rz - rz + gss*lz*rz + hss*nz*rz)*(gss*sz - sz + gss*lz*sz + hss*nz*sz)^2)/(2*(gss + hss - pss*(gss + gu*su*u1u) + su*u1u - qss*(hss + hu*su*u1u) + 1)^(gamma + 3)) - (beta*gamma*vu_u*(Prc - 1)*(gamma + 1)*(gamma + 2)*(gu*pss*su - su + hu*qss*su)^2*(gss*rz - rz + gss*lz*rz + hss*nz*rz))/(2*(gss + hss - pss*(gss + gu*su*u1u) + su*u1u - qss*(hss + hu*su*u1u) + 1)^(gamma + 3)) + (Prc*beta*gamma*su^2*vu_c*(gamma + 1)*(gamma + 2)*(gss*rz - rz + gss*lz*rz + hss*nz*rz))/(2*(gss + hss + phi + su*u1c + 1)^(gamma + 3));
U2b =(beta*gamma*vu_u*(Prc - 1)*(gamma + 1)*(gamma + 2)*(gu*pss*su - su + hu*qss*su)^2*(gb + hb - qss*(hb^2 + gb*ha) - pss*(ga*gb + gb*hb)))/(2*(gss + hss - pss*(gss + gu*su*u1u) + su*u1u - qss*(hss + hu*su*u1u) + 1)^(gamma + 3)) - beta*((Prc*gamma*(gb + hb))/(gss + hss + phi + su*u1c + 1)^(gamma + 1) - (gamma*(Prc - 1)*(gb + hb - qss*(hb^2 + gb*ha) - pss*(ga*gb + gb*hb)))/(gss + hss - pss*(gss + gu*su*u1u) + su*u1u - qss*(hss + hu*su*u1u) + 1)^(gamma + 1)) - beta*sz*((Prc*gamma*gb*sz)/(gss + hss + phi + su*u1c + 1)^(gamma + 1) - (gamma*gb*sz*(Prc - 1))/(gss + hss - pss*(gss + gu*su*u1u) + su*u1u - qss*(hss + hu*su*u1u) + 1)^(gamma + 1) - (Prc*gamma*(gb + hb)*(gamma + 1)*(gss*sz - sz + gss*lz*sz + hss*nz*sz))/(gss + hss + phi + su*u1c + 1)^(gamma + 2) + (gamma*(Prc - 1)*(gamma + 1)*(gss*sz - sz + gss*lz*sz + hss*nz*sz)*(gb + hb - qss*(hb^2 + gb*ha) - pss*(ga*gb + gb*hb)))/(gss + hss - pss*(gss + gu*su*u1u) + su*u1u - qss*(hss + hu*su*u1u) + 1)^(gamma + 2)) - (gamma*pss*(gb*pss + hb*qss - 1))/(gss + hss - gss*pss - hss*qss + 1)^(gamma + 1) - (beta*((Prc*gamma*(gb + hb)*(gamma + 1)*(gamma + 2)*(gss*sz - sz + gss*lz*sz + hss*nz*sz)^2)/(gss + hss + phi + su*u1c + 1)^(gamma + 3) + (2*gamma*gb*sz*(Prc - 1)*(gamma + 1)*(gss*sz - sz + gss*lz*sz + hss*nz*sz))/(gss + hss - pss*(gss + gu*su*u1u) + su*u1u - qss*(hss + hu*su*u1u) + 1)^(gamma + 2) - (gamma*(Prc - 1)*(gamma + 1)*(gamma + 2)*(gss*sz - sz + gss*lz*sz + hss*nz*sz)^2*(gb + hb - qss*(hb^2 + gb*ha) - pss*(ga*gb + gb*hb)))/(gss + hss - pss*(gss + gu*su*u1u) + su*u1u - qss*(hss + hu*su*u1u) + 1)^(gamma + 3) - (2*Prc*gamma*gb*sz*(gamma + 1)*(gss*sz - sz + gss*lz*sz + hss*nz*sz))/(gss + hss + phi + su*u1c + 1)^(gamma + 2)))/2 - (Prc*beta*gamma*su^2*vu_c*(gb + hb)*(gamma + 1)*(gamma + 2))/(2*(gss + hss + phi + su*u1c + 1)^(gamma + 3));
U2a =(beta*gamma*vu_u*(Prc - 1)*(gamma + 1)*(gamma + 2)*(gu*pss*su - su + hu*qss*su)^2*(ga + ha - pss*(ga^2 + gb*ha) - qss*(ga*ha + ha*hb)))/(2*(gss + hss - pss*(gss + gu*su*u1u) + su*u1u - qss*(hss + hu*su*u1u) + 1)^(gamma + 3)) - beta*((Prc*gamma*(ga + ha))/(gss + hss + phi + su*u1c + 1)^(gamma + 1) - (gamma*(Prc - 1)*(ga + ha - pss*(ga^2 + gb*ha) - qss*(ga*ha + ha*hb)))/(gss + hss - pss*(gss + gu*su*u1u) + su*u1u - qss*(hss + hu*su*u1u) + 1)^(gamma + 1)) - beta*sz*((Prc*ga*gamma*sz)/(gss + hss + phi + su*u1c + 1)^(gamma + 1) - (ga*gamma*sz*(Prc - 1))/(gss + hss - pss*(gss + gu*su*u1u) + su*u1u - qss*(hss + hu*su*u1u) + 1)^(gamma + 1) - (Prc*gamma*(ga + ha)*(gamma + 1)*(gss*sz - sz + gss*lz*sz + hss*nz*sz))/(gss + hss + phi + su*u1c + 1)^(gamma + 2) + (gamma*(Prc - 1)*(gamma + 1)*(gss*sz - sz + gss*lz*sz + hss*nz*sz)*(ga + ha - pss*(ga^2 + gb*ha) - qss*(ga*ha + ha*hb)))/(gss + hss - pss*(gss + gu*su*u1u) + su*u1u - qss*(hss + hu*su*u1u) + 1)^(gamma + 2)) - (gamma*pss*(ga*pss + ha*qss - 1))/(gss + hss - gss*pss - hss*qss + 1)^(gamma + 1) - (beta*((Prc*gamma*(ga + ha)*(gamma + 1)*(gamma + 2)*(gss*sz - sz + gss*lz*sz + hss*nz*sz)^2)/(gss + hss + phi + su*u1c + 1)^(gamma + 3) + (2*ga*gamma*sz*(Prc - 1)*(gamma + 1)*(gss*sz - sz + gss*lz*sz + hss*nz*sz))/(gss + hss - pss*(gss + gu*su*u1u) + su*u1u - qss*(hss + hu*su*u1u) + 1)^(gamma + 2) - (gamma*(Prc - 1)*(gamma + 1)*(gamma + 2)*(gss*sz - sz + gss*lz*sz + hss*nz*sz)^2*(ga + ha - pss*(ga^2 + gb*ha) - qss*(ga*ha + ha*hb)))/(gss + hss - pss*(gss + gu*su*u1u) + su*u1u - qss*(hss + hu*su*u1u) + 1)^(gamma + 3) - (2*Prc*ga*gamma*sz*(gamma + 1)*(gss*sz - sz + gss*lz*sz + hss*nz*sz))/(gss + hss + phi + su*u1c + 1)^(gamma + 2)))/2 - (Prc*beta*gamma*su^2*vu_c*(ga + ha)*(gamma + 1)*(gamma + 2))/(2*(gss + hss + phi + su*u1c + 1)^(gamma + 3));
U2ut =beta*((gamma*(Prc - 1)*(gu + hu - pss*(ga*gu + gb*hu) - qss*(gu*ha + hb*hu)))/(gss + hss - pss*(gss + gu*su*u1u) + su*u1u - qss*(hss + hu*su*u1u) + 1)^(gamma + 1) - (Prc*gamma*(gu + hu))/(gss + hss + phi + su*u1c + 1)^(gamma + 1)) - (beta*((Prc*gamma*(gu + hu)*(gamma + 1)*(gamma + 2)*(gss*sz - sz + gss*lz*sz + hss*nz*sz)^2)/(gss + hss + phi + su*u1c + 1)^(gamma + 3) + (2*gamma*gu*sz*(Prc - 1)*(gamma + 1)*(gss*sz - sz + gss*lz*sz + hss*nz*sz))/(gss + hss - pss*(gss + gu*su*u1u) + su*u1u - qss*(hss + hu*su*u1u) + 1)^(gamma + 2) - (2*Prc*gamma*gu*sz*(gamma + 1)*(gss*sz - sz + gss*lz*sz + hss*nz*sz))/(gss + hss + phi + su*u1c + 1)^(gamma + 2) - (gamma*(Prc - 1)*(gamma + 1)*(gamma + 2)*(gu + hu - pss*(ga*gu + gb*hu) - qss*(gu*ha + hb*hu))*(gss*sz - sz + gss*lz*sz + hss*nz*sz)^2)/(gss + hss - pss*(gss + gu*su*u1u) + su*u1u - qss*(hss + hu*su*u1u) + 1)^(gamma + 3)))/2 - beta*sz*((Prc*gamma*gu*sz)/(gss + hss + phi + su*u1c + 1)^(gamma + 1) - (gamma*gu*sz*(Prc - 1))/(gss + hss - pss*(gss + gu*su*u1u) + su*u1u - qss*(hss + hu*su*u1u) + 1)^(gamma + 1) - (Prc*gamma*(gu + hu)*(gamma + 1)*(gss*sz - sz + gss*lz*sz + hss*nz*sz))/(gss + hss + phi + su*u1c + 1)^(gamma + 2) + (gamma*(Prc - 1)*(gamma + 1)*(gu + hu - pss*(ga*gu + gb*hu) - qss*(gu*ha + hb*hu))*(gss*sz - sz + gss*lz*sz + hss*nz*sz))/(gss + hss - pss*(gss + gu*su*u1u) + su*u1u - qss*(hss + hu*su*u1u) + 1)^(gamma + 2)) - (gamma*pss*(gu*pss + hu*qss - 1))/(gss + hss - gss*pss - hss*qss + 1)^(gamma + 1) + (beta*gamma*vu_u*(Prc - 1)*(gamma + 1)*(gamma + 2)*(gu*pss*su - su + hu*qss*su)^2*(gu + hu - pss*(ga*gu + gb*hu) - qss*(gu*ha + hb*hu)))/(2*(gss + hss - pss*(gss + gu*su*u1u) + su*u1u - qss*(hss + hu*su*u1u) + 1)^(gamma + 3)) - (Prc*beta*gamma*su^2*vu_c*(gu + hu)*(gamma + 1)*(gamma + 2))/(2*(gss + hss + phi + su*u1c + 1)^(gamma + 3));
U2zt =beta*((Prc*gamma*(gss*rz - rz + gss*lz*rz + hss*nz*rz))/(gss + hss + phi + su*u1c + 1)^(gamma + 1) - (gamma*(Prc - 1)*(gss*rz - rz + gss*lz*rz + hss*nz*rz))/(gss + hss - pss*(gss + gu*su*u1u) + su*u1u - qss*(hss + hu*su*u1u) + 1)^(gamma + 1)) - (vu_c*((Prc*beta*gamma*rz*su^2*(gamma + 1))/(gss + hss + phi + su*u1c + 1)^(gamma + 2) - (Prc*beta*gamma*su^2*(gamma + 1)*(gamma + 2)*(gss*rz - rz + gss*lz*rz + hss*nz*rz))/(gss + hss + phi + su*u1c + 1)^(gamma + 3)))/2 + (vu_u*((beta*gamma*rz*(Prc - 1)*(gamma + 1)*(gu*pss*su - su + hu*qss*su)^2)/(gss + hss - pss*(gss + gu*su*u1u) + su*u1u - qss*(hss + hu*su*u1u) + 1)^(gamma + 2) - (beta*gamma*(Prc - 1)*(gamma + 1)*(gamma + 2)*(gu*pss*su - su + hu*qss*su)^2*(gss*rz - rz + gss*lz*rz + hss*nz*rz))/(gss + hss - pss*(gss + gu*su*u1u) + su*u1u - qss*(hss + hu*su*u1u) + 1)^(gamma + 3)))/2 - (beta*((gamma*(Prc - 1)*(gamma + 1)*(gamma + 2)*(gss*rz - rz + gss*lz*rz + hss*nz*rz)*(gss*sz - sz + gss*lz*sz + hss*nz*sz)^2)/(gss + hss - pss*(gss + gu*su*u1u) + su*u1u - qss*(hss + hu*su*u1u) + 1)^(gamma + 3) - (Prc*gamma*(gamma + 1)*(gamma + 2)*(gss*rz - rz + gss*lz*rz + hss*nz*rz)*(gss*sz - sz + gss*lz*sz + hss*nz*sz)^2)/(gss + hss + phi + su*u1c + 1)^(gamma + 3)))/2 - lz/(gss + hss - gss*pss - hss*qss + 1)^gamma + beta*rz*((Prc - 1)/(gss + hss - pss*(gss + gu*su*u1u) + su*u1u - qss*(hss + hu*su*u1u) + 1)^gamma - Prc/(gss + hss + phi + su*u1c + 1)^gamma) + beta*sz*((gamma*(Prc - 1)*(gamma + 1)*(gss*rz - rz + gss*lz*rz + hss*nz*rz)*(gss*sz - sz + gss*lz*sz + hss*nz*sz))/(gss + hss - pss*(gss + gu*su*u1u) + su*u1u - qss*(hss + hu*su*u1u) + 1)^(gamma + 2) - (Prc*gamma*(gamma + 1)*(gss*rz - rz + gss*lz*rz + hss*nz*rz)*(gss*sz - sz + gss*lz*sz + hss*nz*sz))/(gss + hss + phi + su*u1c + 1)^(gamma + 2)) - (beta*rz*((Prc*gamma*(gamma + 1)*(gss*sz - sz + gss*lz*sz + hss*nz*sz)^2)/(gss + hss + phi + su*u1c + 1)^(gamma + 2) - (gamma*(Prc - 1)*(gamma + 1)*(gss*sz - sz + gss*lz*sz + hss*nz*sz)^2)/(gss + hss - pss*(gss + gu*su*u1u) + su*u1u - qss*(hss + hu*su*u1u) + 1)^(gamma + 2)))/2 - (gamma*pss*(gss + gss*lz + hss*nz - 1))/(gss + hss - gss*pss - hss*qss + 1)^(gamma + 1);

%Collect the equations
y = [ZOb;ZOa;U1;U2;U1b;U1a;U1ut;U1zt;U2b;U2a;U2ut;U2zt];