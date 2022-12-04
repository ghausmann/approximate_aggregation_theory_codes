%The following function computes the stationary distribution, iterating
%over the asset transition matrix. The credit goes to Guerierri and
%Lorenzoni, from the replication codes of their 2017 QJE article. 
%The inputs are b_grid (the asset grid with length n_b), Pr (Probability
%transition matrix with size n_z-by-n_z), b_pold (the savings rule,a matrix
%with size n_z-by-n_b), D0 (an initial guess for the distribution, a
%n_z-by-n_b matrix), and a tolerance value.

function y = ss_distribution(b_grid,Pr,b_pold,D0,mytol)

msize = size(b_pold);

n_z = msize(1);
n_b = msize(2);

% Young's lottery method of interpolation
[~,~ ,ib_pol] = histcounts(b_pold, b_grid);
wei = (b_pold - b_grid(ib_pol)) ./ (b_grid(ib_pol+1) - b_grid(ib_pol));

dif = 1;
D  = D0;

while dif > mytol
    Di = zeros(n_z, n_b);
    
    %Reallocate mass of households from a given distribution according to the optimal saving rule 
    for s = 1:n_z
        for i = 1:n_b

            Di(s, ib_pol(s, i))     = (1 - wei(s, i))  * D(s, i) + Di(s, ib_pol(s, i));
            Di(s, ib_pol(s, i) + 1) = wei(s, i)        * D(s, i) + Di(s, ib_pol(s, i) + 1);

        end
    end
    % make sure that distribution integrates to 1
    Di = Di / sum(sum(Di));
    Di = (Pr')*Di; %account for future individual luck
    % check convergence
    dif = max(max(abs(Di - D)));
    D = Di;
end

y = D;