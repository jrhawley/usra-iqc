%---Effect Discrete Wigner function
% Inputs:
%   E = n-by-n-by-n array of measurement projectors
% Outputs:
%   W = discrete wigner representation
%   M = mana

function [W, M] = ewigner(E)
    [~,~,n] = size(E);
    global A
    W = zeros(n,n,n);
    for j = 1:n
        for k = 1:n
            for l = 1:n
                W(j,k,l) = trace(E(:,:,l)*A(:,:,j,k)); %W_E_l = Tr(E_l*A_jk)
            end
        end
    end
    
    M = zeros(n,n);
    for j = 1:n
        for k = 1:n
            for l = 1:n
                M(j,k) = M(j,k) + abs(W(j,k,l));
            end
        end
    end
end