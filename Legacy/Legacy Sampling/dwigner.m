%---Density Discrete Wigner function
% Inputs:
%   rho = n-by-n density operator
% Outputs:
%   W = discrete wigner representation
%   M = mana

function [W, M] = dwigner(rho)
    [~,n] = size(rho);
    global A
    W = zeros(n);
    for j = 1:n
        for k = 1:n
            W(j,k) = 1/n*trace(rho*A(:,:,j,k));
        end
    end
    
    M = sum(sum(abs(W))); % return mana of state matrix
end