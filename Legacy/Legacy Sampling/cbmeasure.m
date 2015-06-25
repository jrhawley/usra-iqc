function o = cbmeasure (rho)
    [n,~] = size(rho);
    p = zeros(n,1);
    
    % build probabilities
    for j = 1:n
        E = zeros(n);
        E(j,j) = 1;
        p(j) = trace(E*rho);
    end
    
    % randomly select measurement
    o = randsample(n, 1, true, p);
    o = mod(o - 1, n); % 0-index outcomes
end