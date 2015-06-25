%---indlm1
% Converts linear index to indices of an n-by-n grid (1-indexed)
% Inputs:
%   n = dimension of square grid
%   k = row or column vector, with values from 1 to n^2
% Outputs:
%   v = matrix where the ith row contains indices of ith entry of k
function v = indlm1(n,k)
    %  determine column or row vector and turn k into column vector
    [rk,ck] = size(k);
    if 1 == rk
        nk = ck;
        k = k';
    elseif 1 == ck
        nk = rk;
    else
        error('Bad dimensions on k');
    end
    
    % make sure all values are within 1 to n^2
    if sum(n^2 < k) > 0
        error('Illegal value > n^2 in k');
    elseif sum(0 >= k) > 0
        error('Illegal value < 1 in k');
    end
    
    v = zeros(nk,2);
    v(:,2) = mod(k,n);
    for j = 1:nk
        if 0 == v(j,2)
            v(j,2) = n;
        end
    end
    v(:,1) = (k - v(:,2))/n + 1;
end
    