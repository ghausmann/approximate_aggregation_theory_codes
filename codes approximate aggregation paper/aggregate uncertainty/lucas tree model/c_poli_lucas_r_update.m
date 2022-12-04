%This function uses the endogenous grid method, mapping an old consumption
%rule to a new one, accounting for both individual and aggregate risk. It
%assumes that the equilibrium asset price is a linear polynomial of the
%type p = pss + lY*(Y-1), where Y output, and (pss,ly) two coefficients.
%Inputs are P (vector of parameters), L (the vector [pss,lY]) a_grid and zi
%(grids for asset of length n_a and productivity values of length n_z), Pr
%(Probability transition matrix of size n_z-by_n_z), Yg (the grid for
%output of size n_Y), Pr_Y (Probability transition matrix of size
%n_Y-by-n_Y), c_pol (the old/tomorrow's consumption rule, a
%n_z-by-n_a-by-n_Y array), and order (the other of approximation).

function y = c_poli_lucas_r_update(P,L,a_grid,zi,Pr,Yg,PrY,c_pol,order)

%Model parameters
beta = P(1);
gamma = P(2);
alpha = P(3);

pss = L(1);
lY = L(2);
%Include the quadratic term except if order==1
if order ==1
    lY2 = 0;
else
    lY2 = L(3);
end

%lenghts of the asset, productivity, and output grids
n_a = length(a_grid);
n_z = length(zi);
n_Y = length(Yg);

c_poli = zeros(n_z,n_a,n_Y);

%Compute each of today's conditional consumption rules (SxI matrix)
for iy = 1:n_Y

    %select the proper conditional distribution for tomorrow's output
    Prz_cond = PrY(iy,:);
    Prz_cond_res = reshape(Prz_cond,1,1,n_Y);

    %array for conditional expected marginal utilities
    Prc_pol = zeros(n_z,n_a,n_Y);

    %Compute each expected marginal utility, conditional on z
    for iiy=1:n_Y

        Prc_pol(:,:,iiy) = ((pss + lY*(Yg(iiy)-1) + 0.5*lY2*(Yg(iiy)-1)^2 ) + alpha*( Yg(iiy) ))*Pr*(1./c_pol(:,:,iiy)).^gamma;

    end

    %Compute the unconditional expected marginal utility
    sPrc_cpol = sum(Prz_cond_res.*Prc_pol,3);

    %Compute today's consumption using the Euler equation, for fixed values of
    %today's savings and the asset price (p = p_ss + ly*z).
    c_i = ((1/(  (beta/(pss + lY*(Yg(iy)-1) + 0.5*lY2*(Yg(iy)-1)^2 )  )))./(     sPrc_cpol   )     ).^(1/gamma);

    %Compute today's implied initial level of assets.
    a_grid1 = ( c_i + (pss + lY*(Yg(iy)-1)  + 0.5*lY2*(Yg(iy)-1)^2 )*(1+a_grid) - (1+zi')*(1-alpha)*(Yg(iy)) )/((pss + lY*(Yg(iy)-1) + 0.5*lY2*(Yg(iy)-1)^2 ) + alpha*(Yg(iy))) - 1;

    c_poliy = zeros(n_z,n_a);

    %Interpolate to complete the mapping from the asset grid to the new
    %consumption rule
    for j=1:n_z

        c_poliy(j,:) = interp1(a_grid1(j,:), c_i(j,:), a_grid, 'linear', 'extrap');

    end

    %Correct values that violate the constraints (-1<=a_poli<=a_grid(I)):

    %Compute the savings rule
    a_poliy = (  (1+zi')*(1-alpha)*(Yg(iy)) + ((pss + lY*(Yg(iy)-1) + 0.5*lY2*(Yg(iy)-1)^2 )+alpha*(Yg(iy)))*(1+a_grid)   -  c_poliy )/(pss + lY*(Yg(iy)-1) + 0.5*lY2*(Yg(iy)-1)^2 ) - 1;
    %Consumption if the borrowing constraint binds
    c_poli_d = (1+zi')*(1-alpha)*(Yg(iy)) + ((pss + lY*(Yg(iy)-1) + 0.5*lY2*(Yg(iy)-1)^2 )+alpha*(Yg(iy)))*(1+a_grid);
    %Consumption if the upper limit for savings binds
    c_poli_u = (1+zi')*(1-alpha)*(Yg(iy)) + ((pss + lY*(Yg(iy)-1) + 0.5*lY2*(Yg(iy)-1)^2 )+alpha*(Yg(iy)))*(1+a_grid) - (pss + lY*(Yg(iy)-1) + 0.5*lY2*(Yg(iy)-1)^2 )*(1+a_grid(n_a));

    %Replace values that violate constraints by the proper ones
    c_poliy(a_poliy<=-1) = c_poli_d(a_poliy<=-1);
    c_poliy(a_poliy>=a_grid(n_a)) = c_poli_u(a_poliy>=a_grid(n_a));

    %assign the resulting conditional consumption rule (n_z-by-n_a) in the new
    %consumption array (n_z-by-n_a-by-n_Y)
    c_poli(:,:,iy) = c_poliy;

end

%Return today's policy rule
y = c_poli;