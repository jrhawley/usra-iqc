%---indml0
% Converts indices of an n-by-n grid to linear index (0-indexed)
% Inputs:
%   n = dimension of square grid
%   v = nk-by-2 matrix of coordinates
% Outputs:
%   k = nk-by-1 vector of linear indices
function k = indml0(n,v)
    if sum(sum(v > n-1)) > 0
        error('Illegal value > n-1');
    elseif sum(sum(v < 0)) > 0
        error('Illegal value < 0');
    end
    
    k = n*(v(:,1)) + v(:,2);
end