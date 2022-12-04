%This function computes the euler errors of the linear approximations
%to the policy and asset rules in the two-asset model.
%The inputs are , P (vector of %parameters), x (the 12 coefficients
%characterizing the solution), and S (current states).
%The output is the unit-free euler error.

function y = local_2assets_errors(P,x,S)

%Paremeters
beta = P(1);
gamma = P(2);
phi = P(3);
su = P(4);
sz = P(5);
rz = P(6);

%Coefficients
qss = x(1);
lz = x(2);
pss = x(3);
nz = x(4);
gb = x(5);
ga = x(6);
gu = x(7);
hb = x(8);
ha = x(9);
hu = x(10);
gss = x(11);
hss = x(12);

%coefficients of the discriminant function
fss = qss*hss + pss*gss;
fb = qss*hb + pss*gb;
fa = qss*ha + pss*ga;
fu = qss*hu + pss*gu;

%states
a_1 = S(1);
b_1 = S(2);
ut = S(3);
zt = S(4);

%asset prices today
qt = qss + lz*zt;
pt = pss + nz*zt;

%saving rules today
bt = hss + ha*(a_1-gss)  + hb*(b_1-hss)  + hu*ut;
at = gss + ga*(a_1-gss)  + gb*(b_1-hss)  + gu*ut;

%consumption today
ct = 1 + ut + zt + (1 - zt)*a_1 + b_1 - ( pss*at + gss*(pt-pss) ) - ( qss*bt + hss*(qt-qss) ) ;

%saving rules tomorrow
ggt1 = gss + ga*(at-gss)  + gb*(bt-hss);
hht1 = hss + ha*(at-gss)  + hb*(bt-hss);
fft1 = fss + fa*(at-gss)  + fb*(bt-hss);

at1u = @(u1)ggt1 + gu*su*u1;
bt1u = @(u1)hht1 + hu*su*u1;

%exogenous processes tomorrow
ut1 = @(u1)su*u1;
zt1 = @(v1)rz*zt + sz*v1;

%asset prices tomorrow
qbt1 = @(v1)qss + lz*zt1(v1);
qat1 = @(v1)pss + nz*zt1(v1);

%consumption tomorrow, unconstrained and contrained
ct1u = @(u1,v1)1 + ut1(u1) + zt1(v1) + (1 - zt1(v1))*at + bt - ( pss*at1u(u1) + gss*(qat1(v1)-pss) ) - ( qss*bt1u(u1) + hss*(qbt1(v1)-qss) ) ;
ct1c = @(u1,v1)1 + ut1(u1) + zt1(v1) + (1 - zt1(v1))*at + bt - (  gss*(qat1(v1)-pss) ) - (  hss*(qbt1(v1)-qss) ) + phi;

%Bivariate normal distribution
npdf = @(u1,v1)(((0.5/pi))*exp(-0.5*(u1.^2 + v1.^2)));

%Marginal utilities tomorrow
%For the Euler equation of the risky asset, constrained and unconstrained
%cases
myexp_int_a_U1 = @(u1,v1)( (1 - zt1(v1)).*(ct1c(u1,v1).^-gamma)  ).*npdf(u1,v1);
myexp_int_b_U1 = @(u1,v1)( (1 - zt1(v1)).*(ct1u(u1,v1).^-gamma)  ).*npdf(u1,v1);
%For the Euler equation of the safe asset, constrained and unconstrained
%cases
myexp_int_a_U2 = @(u1,v1)( (ct1c(u1,v1).^-gamma)  ).*npdf(u1,v1);
myexp_int_b_U2 = @(u1,v1)( (ct1u(u1,v1).^-gamma)  ).*npdf(u1,v1);

%cutoff point (u1<mycut_1 is the constrained case)
mycut_1 = - (fft1 + phi )./(fu*su);
mycut1 = mycut_1 + 1e-5;

%Integrate to compute the expectations of
myexp_a_U1 = integral2(myexp_int_a_U1,-6,mycut_1,-6,6);
myexp_b_U1 = integral2(myexp_int_b_U1,mycut1,6,-6,6);
myexp_a_U2 = integral2(myexp_int_a_U2,-6,mycut_1,-6,6);
myexp_b_U2 = quad2d(myexp_int_b_U2,mycut1,6,-6,6);

%overall expected marginal utilities
myexp_U1 = myexp_a_U1 + myexp_b_U1;
myexp_U2 = myexp_a_U2 + myexp_b_U2;

%compute the unit-free euler errors
errors_U1 = abs( ((beta*myexp_U1/pt)^(-1/gamma))/ct -1  );
errors_U2 = abs( ((beta*myexp_U2/qt)^(-1/gamma))/ct -1  );

y = [errors_U1;errors_U2];
