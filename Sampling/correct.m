%---Correction Function
% turns values close to -1, 0, 1 (ie within a small range) to exactly those
% numbers
% Inputs:
%   x = value
% Outputs:
%   X = corrected value
function X = correct(x)
    if (abs(x-1) < 10*eps)
        X = 1;
    elseif (abs(x-0) < 10*eps)
        X = 0;
    elseif (abs(x+1) < 10*eps)
        X = -1;
    else
        X = x;
    end
end