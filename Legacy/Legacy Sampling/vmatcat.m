function R = vmatcat(M)
    [~,n] = size(M);
    R = zeros(n^2,1);
    for j = 0:n-1
        R((j*n+1):(j+1)*n) = M(:,j+1);
    end
end