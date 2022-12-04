%This function calculates the M(1) approximation.
%As in the RA case, it takes the form
%(q-q_ss) = lZ*(Z-Zss)
%where q_ss is the SS bond price, Z a exogenous process, and lZ
%the marginal response of the bond price to changes in Z.
%In this case we have a explicit theoretical formula for lZ, which uses as
%inputs the sequences of average derivatives. 
%The inputs are rho (persistence of the exogenous process), Tmax (maximum
%number of periods), idgdq (average partial derivatives w.r.t. future bond
%prices), and idgdD (average partial derivatives w.r.t. future exogenous
%shocks).

function y = savings_M1_Z(rho,Tmax,idgdq,idgdZ)


mytime = 0:1:(Tmax-1);

%The formula comes from the equilibrium condition dBp'dZ = 0, where dBp'dZ
%is the average (total) derivative of today's savings rule w.r.t today's
%exogenous shock. The condition says that today's market clearing must hold
%despite a small shock to the exogenous process). Here dBp'dZ depends on lZ,
%and it can be decomposed using the average partial derivatives, from which
%we can solve for lZ. 
mynum = sum(  idgdZ.*(rho.^mytime)   );
myden = sum(  idgdq.*(rho.^mytime)   );

y = -mynum/myden;