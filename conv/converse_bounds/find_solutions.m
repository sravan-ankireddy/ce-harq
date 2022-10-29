function [nsol, Sols] = find_solutions(n)

% The function aims at identifying all non-negative integer solutions to
% the Diophantine equation
%   k_1 + 2*k_2 + ... +n*k_n = n.
% We follow the algorithm described in [1] to solve this problem.
%
% References:
% [1] S. Blinnikov and R. Moessner, "Expansions for nearly Gaussian 
%     distributions", Astron. Astrophys. Suppl. Ser., 1998.
%
% Input parameters:
%   1) n: a positive integer
%
% Output parameters:
%   1) nsol: a scalar denoting the # solutions
%   2) Sols: a matrix with each row denoting a solution {k_1, k_2, ..., k_n}
%
%
% Written by Hengjie Yang (hengjie.yang@ucla.edu)   01/11/22
%


% basic parameters
% r = 3; % the base used to represent k_n*r^(n-1) + ... + k_2*r + k_1
ks = zeros(1, n); % the ordering in 'ks' is [k_1, k_2, ..., k_n]
ks(1) = n;

% basis = zeros(1, n); 
% for ii = 1:n
%     basis(ii) = r^(n-ii);
% end
mold = 1;
nsol = 1;
Sols = [];
Sols = [Sols; ks];


while mold < n
    m = 1;
    sumcur = n;
    while 1
        sumcur = sumcur - ks(m)*m + (m + 1);
        ks(m) = 0;
        ks(m + 1) = ks(m + 1) + 1;
        m = m + 1;
        if (sumcur <= n) || (m > mold)
            break
        end
    end
    if m > mold
        mold = m;
    end
    ks(1) = n - sumcur;
    nsol = nsol + 1;
    Sols = [Sols; ks];
end











end