%This function calculates lZ, from the implied law of motion Kp = K_ss +
%lK1*(K-Kss) + lZ*(Z-1) assumed by households, using the explicit formula
%from the paper.
%Inputs are P (vector of coefficients), r (SS rental price), T (maximum
%number of periods), idgdX (average partial derivatives w.r.t. future
%aggregate inputs), and lK (already known).

function y = ks_M1_Z(P,r,Tmax,idgdX,lK)

%Parameter values
alpha = P(3);
rho = P(5);

K = (alpha/r)^(1/(1-alpha)); %SS aggregate capital
w = (1-alpha)*(K^alpha)*((1-alpha)); %SS wage rate
%Derivatives of rental price w.r.t. current states, evaluated at SS
rK = -(1-alpha)*r/K;
wK = alpha*w/K;
rZ =r;
wZ =w;

%Sequences of average derivatives
idkdr = idgdX(1,:); %w.r.t future rental prices
idkdw = idgdX(2,:); %w.r.t future wage rates
idkdK = idgdX(3,:); %w.r.t. the wage rate

%Compute lZ using the corresponding formula (see the paper for details).
S_my_product = zeros(1,Tmax);

for t=0:(Tmax-1)
   
    mytime = 0:1:t;
    S_my_product(t+1) = sum( (rho.^(t-(mytime))).*(lK.^(mytime)) );
        
end

Kt_Z = [0 S_my_product(1:Tmax-1)];

nmytime=0:1:(Tmax-1);
my_num = -sum( (rZ*idkdr + wZ*idkdw).*(rho.^nmytime) );
my_den = sum(  ( rK*idkdr + wK*idkdw + idkdK).*Kt_Z    );

%Return lZ
y = my_num/my_den;

