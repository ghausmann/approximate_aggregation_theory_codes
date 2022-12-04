%This function computes the average partial derivatives (2nd and 3rd order)
%of functions of today's savings rule with respect to current individual
%assets.
%Inputs are b_grid (asset grid with length n_b), b_poli (SS savings rule
%with size n_u-by-n_b), Pr (the probability transition matrix with size
%n_u-by-n_u), and Dss (stationary wealth distribution with size
%n_u-by-n_b).

function y = savings_idgds(b_grid,b_poli,Pr,Dss)

mDss = sum(Dss); %Marginal wealth distribution

%Functions of the savings rule, as expected values over ui, conditional on b: hp_n = E(bp^n|b)
hp_1 = Pr*b_poli; %For the mean
hp_2 = Pr*(b_poli.^2); %For the variance
hp_3 = Pr*(b_poli.^3); %For the "skewness"

%Because ui is purely transitory, all the rows are indentical. So just pick
%the first one.
hp_1 = hp_1(1,:);
hp_2 = hp_2(1,:);
hp_3 = hp_3(1,:);

%Numerical derivatives w.r.t. current asset holdings
%For E(bp)
dhp_1 = gradient(hp_1,b_grid); %first derivative
d2hp_1 = gradient(dhp_1,b_grid); %second derivative
d3hp_1 = gradient(d2hp_1,b_grid); %third derivative
%For E(bp^2)
dhp_2 = gradient(hp_2,b_grid);
d2hp_2 = gradient(dhp_2,b_grid);
d3hp_2 = gradient(d2hp_2,b_grid);
%For E(bp^3)
dhp_3 = gradient(hp_3,b_grid);
d2hp_3 = gradient(dhp_3,b_grid);
d3ph_3 = gradient(d2hp_3,b_grid);

%Compute the average derivatives, properly weighted.
dBp_dV_ind = (0.5)*(  sum( mDss.*d2hp_1    )   );
dVp_dV_ind = (0.5)*( sum(mDss.*d2hp_2  )   );
dM3p_dV_ind = (0.5)*(  sum(mDss.*d2hp_3  )   );

dBp_dM3_ind = (1/6)*(  sum( mDss.*d3hp_1    )   );
dVp_dM3_ind = (1/6)*(  sum( mDss.*d3hp_2    )   );
dM3p_dM3_ind = (1/6)*(  sum( mDss.*d3ph_3    )   );


y = [dBp_dV_ind dVp_dV_ind dM3p_dV_ind;dBp_dM3_ind dVp_dM3_ind dM3p_dM3_ind];