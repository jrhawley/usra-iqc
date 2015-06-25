%---indml1
% Converts indices of an n-by-n grid to linear index (1-indexed)
% Inputs:
%   n = dimension of square grid
%   v = nk-by-2 matrix of coordinates
% Outputs:
%   k = nk-by-1 vector of linear indices
function k = indml1(n,v)
    if sum(sum(v > n)) > 0
        error('Illegal value > n');
    elseif sum(sum(v < 1)) > 0
        error('Illegal value < 1');
    end
    
    k = n*(v(:,1)-1) + v(:,2);
end