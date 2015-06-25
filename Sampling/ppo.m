function A = ppo(n)
    A = zeros(n,n,n,n);
    if (n == 2)
        X = [0 1; 1 0];
        Z = [1 0; 0 -1];
        Y = 1i*X*Z;
        I = eye(2);
        
        for j = 1:n
            for k = 1:n
                A(:,:,j,k) = 1/2*((-1)^(j-1)*Z + (-1)^(k-1)*X + (-1)^(j+k-2)*Y + I);
            end
        end
    else
        for j = 0:n-1
            for k = 0:n-1
                for l = 0:n-1
                    for m = 0:n-1
                        A(l+1,m+1,j+1,k+1) = kdelta(mod(2*j, n), mod(l+m, n))*exp(2*pi*1i/n*k*(l-m));
                    end
                end
            end
        end
    end
end