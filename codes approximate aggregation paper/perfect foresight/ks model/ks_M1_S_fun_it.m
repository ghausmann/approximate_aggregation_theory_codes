%This function takes as an input a given value lK1 from the implied law of
%motion Kp = K_ss + lK1*(K-Kss) assumed by households, and
%updates it using time iteration, thus returning a new value lK0.
%Inputs are P (vector of coefficients), r (SS rental price), T (maximum
%number of periods), idgdK (average derivatives w.r.t. future prices), h
%(small shock to today's capital), and lK1.

function y = ks_M1_S_fun_it(P,r,T,idgdX,h,lK1)

%Model parameters
alpha = P(3);

K = (alpha/r)^(1/(1-alpha)); %SS aggregate capital
w = (1-alpha)*(K^alpha)*((1-alpha)); %SS wage rate
%SS partial derivatives of prices w.r.t. aggregate capital
rK = -(1-alpha)*r/K;
wK = alpha*w/K;

idgdr = idgdX(1,:); %w.r.t the rental price
idgdw = idgdX(2,:); %w.r.t. the wage rate
idgdX = idgdX(3,:); %w.r.t. the wage rate
%--------------------------------------------------------------------------
%K shock
%--------------------------------------------------------------------------
Kt = zeros(1,T);
Kt(1) = h; %Small shock to today's aggregate capital

for t = 2:T

    Kt(t) =lK1*(Kt(t-1));

end
%Derivative of the whole sequence w.r.t. today's shock (in this simple
%case, we know it is the sequence {1,lK1,lK1^2,...}).
dKdK_1 = Kt/h;

%Construct the implied lK0, using the formula from the paper.
%The key time-iteration idea: just for today dKt/dKt_1=lK0, and for any
%future period, dKt/dKt_1=lK1
md = sum( (idgdr(2:T)*rK+ idgdw(2:T)*wK + idgdX(2:T)).*(dKdK_1(1:T-1)) );
lK0 = -(idgdr(1)*rK+ idgdw(1)*wK + idgdX(1))/md;

y = lK0;