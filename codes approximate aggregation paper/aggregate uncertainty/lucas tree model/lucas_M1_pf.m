%This function calculates the M(1) approximation. As in the RA case, it
%takes the form p = p_ss + lY*(Y-1), where p_ss is the SS asset price, Y
%current output, and lY the marginal response of the asset price to changes
%in Y. In this case we have a explicit theoretical formula for lY, which
%uses as inputs the sequences of average derivatives.
%The inputs are rho (persistence of the output process), T (maximum number
%of periods), and idgdX (average partial derivatives w.r.t. future output
%and asset prices).

function y = lucas_M1_pf(rho,T,idgdX)

idgdY = idgdX(1,:);
idgdp = idgdX(2,:);

mytime = 0:1:(T-1);

%The formula comes from the equilibrium condition dAp'dZ = 0, where dAp'dZ
%is the average (total) derivative of today's asset rule w.r.t today's
%output. The condition says that today's market clearing must hold
%despite a small shock to output). Here dAp'dZ depends on lY,
%and it can be decomposed using the average partial derivatives, from which
%we can solve for lY.
mynum = sum(  idgdY.*(rho.^mytime)   );
myden = sum(  idgdp.*(rho.^mytime)   );

y = -mynum/myden;