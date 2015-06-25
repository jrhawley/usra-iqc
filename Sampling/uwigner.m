%---Unitary Discrete Wigner function
% Inputs:
%   U = n-by-n unitary operator
% Outputs:
%   W = discrete wigner representation
%   M = mana

function [W, M] = uwigner(U)
    [~,n] = size(U);
    global A
    W = zeros(n,n,n,n);
    for j = 1:n % (j,k) are components of beta
        for k = 1:n
            for l = 1:n % (l,m) are components of alpha
                for m = 1:n
                    W(j,k,l,m) = 1/n*trace(A(:,:,j,k)*U*A(:,:,l,m)*conj(U'));
                end
            end
        end
    end
    
    M = zeros(n,n);
    for j = 1:n
        for k = 1:n
            for l = 1:n
                for m = 1:n
                    M(l,m) = M(l,m) + abs(W(j,k,l,m));
                end
            end
        end
    end
end