%This script derives 10 out of the 12 equations that solves the two-asset model
%using symbolic differentiation.

clear;clc;

syms beta gamma phi su sz rz; %parameters
syms a_1 b_1 ut zt; %states
syms u1u u1c v1; %future random shocks (u1u=future idiosyncratic shock if unconstrained)
%Coefficients saving rules
syms gss ga gb gu; %for risky asset at
syms hss ha hb hu; %for safe asset bt

syms qss nz pss lz; %Coefficients asset prices
syms Prc vu_u vu_c; %probability objects (Prc=prob of constrained, vu_u=variance of u1u)

%asset prices today
pt = pss + lz*zt;
qt = qss + nz*zt;

%saving rules today
at = gss + ga*(a_1-gss)  + gb*(b_1-hss)  + gu*ut;
bt = hss + ha*(a_1-gss)  + hb*(b_1-hss)  + hu*ut;

%consumption today (It uses the linear approximations pt*at = pss*at +
%gss*(pt-pss) and qt*bt = qss*bt + hss*(qt-qss))
%
ct = 1 + ut + zt + (1 - zt)*a_1 + b_1 - ( pss*at + gss*(pt-pss) ) - ( qss*bt + hss*(qt-qss) ) ;

%exogenous processes tomorrow
zt1 = rz*zt + sz*v1; %future aggregate productivity
ut1u = su*u1u; %conditional on unconstrained
ut1c = su*u1c; %conditional on constrained

%asset prices tomorrow
pt1 = pss + lz*zt1;
qt1 = qss + nz*zt1;

%saving rules tomorrow (if unconstrained)
at1u = gss + ga*(at-gss)  + gb*(bt-hss) + gu*ut1u;
bt1u = hss + ha*(at-gss)  + hb*(bt-hss) + hu*ut1u;

%consumption tomorrow
ct1u = 1 + ut1u + zt1 + (1 - zt1)*at + bt - ( pss*at1u + gss*(pt1-pss) ) - ( qss*bt1u + hss*(qt1-qss) ); %if unconstrained
ct1c = 1 + ut1c + zt1 + (1 - zt1)*at + bt - ( gss*(pt1-pss) ) - ( hss*(qt1-qss) ) + phi; %if constrained

%Residual functions
U1 = Prc*(ct1c^(-gamma))*beta + (1-Prc)*(ct1u^(-gamma))*beta - qt*(ct^(-gamma)); %Euler equation for safe asset
U2 = beta*( Prc*(ct1c^(-gamma)) + (1-Prc)*(ct1u^(-gamma)) )*(1 - zt1) - pt*(ct^(-gamma)); %Euler equation for risky asset

%Derivatives w.r.t. the states
U1a = diff(U1,a_1);
U1b = diff(U1,b_1);
U1ut = diff(U1,ut);
U1zt = diff(U1,zt);

U2a = diff(U2,a_1);
U2b = diff(U2,b_1);
U2ut = diff(U2,ut);
U2zt = diff(U2,zt);

%Evaluate at the steady-state
a_1 = gss;  b_1 = hss; ut=0 ; zt = 0;

U1 = subs(U1);
U1a = subs(U1a);
U1b = subs(U1b);
U1ut = subs(U1ut);
U1zt = subs(U1zt);

U2 = subs(U2);
U2a = subs(U2a);
U2b = subs(U2b);
U2ut = subs(U2ut);
U2zt = subs(U2zt);

%Approximate expectations using the delta method
U1 = U1 + 0.5*(diff(diff(U1,v1),v1) + diff(diff(U1,u1u),u1u)*vu_u + diff(diff(U1,u1c),u1c)*vu_c );
U2 = U2 + 0.5*(diff(diff(U2,v1),v1) + diff(diff(U2,u1u),u1u)*vu_u + diff(diff(U2,u1c),u1c)*vu_c );
U1b = U1b + 0.5*(diff(diff(U1b,v1),v1) + diff(diff(U1b,u1u),u1u)*vu_u + diff(diff(U1b,u1c),u1c)*vu_c );
U1a = U1a + 0.5*(diff(diff(U1a,v1),v1) + diff(diff(U1a,u1u),u1u)*vu_u + diff(diff(U1a,u1c),u1c)*vu_c );
U1ut = U1ut + 0.5*(diff(diff(U1ut,v1),v1) + diff(diff(U1ut,u1u),u1u)*vu_u + diff(diff(U1ut,u1c),u1c)*vu_c );
U1zt = U1zt + 0.5*(diff(diff(U1zt,v1),v1) + diff(diff(U1zt,u1u),u1u)*vu_u + diff(diff(U1zt,u1c),u1c)*vu_c );

U2b = U2b + 0.5*(diff(diff(U2b,v1),v1) + diff(diff(U2b,u1u),u1u)*vu_u + diff(diff(U2b,u1c),u1c)*vu_c );
U2a = U2a + 0.5*(diff(diff(U2a,v1),v1) + diff(diff(U2a,u1u),u1u)*vu_u + diff(diff(U2a,u1c),u1c)*vu_c );
U2ut = U2ut + 0.5*(diff(diff(U2ut,v1),v1) + diff(diff(U2ut,u1u),u1u)*vu_u + diff(diff(U2ut,u1c),u1c)*vu_c );
U2zt = U2zt + 0.5*(diff(diff(U2zt,v1),v1) + diff(diff(U2zt,u1u),u1u)*vu_u + diff(diff(U2zt,u1c),u1c)*vu_c );


%Evaluate at the zero mean of aggregate innovation, and return the
%equations
v1 = 0;

U1 = subs(U1)
U2 = subs(U2)

U1b = subs(U1b)
U1a = subs(U1a)
U1ut = subs(U1ut)
U1zt = subs(U1zt)


U2b = subs(U2b)
U2a = subs(U2a)
U2ut = subs(U2ut)
U2zt = subs(U2zt)





