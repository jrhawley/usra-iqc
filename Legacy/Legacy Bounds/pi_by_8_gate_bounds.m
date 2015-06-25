p0 = 50;
p = primes(p0);
[~,np] = size(p);
mana = ones(1,np);

%---d = 2------------------------------------------------------------------
U = diag([1,exp(1i*pi/4)]);
Ut = U'; % ' automatically implies conjugate transpose
X = [0 1; 1 0];
Z = [1 0; 0 -1];
Dxz = zeros(2,2,2,2);
Dxzt = zeros(2,2,2,2);
Dxz(:,:,1,1) = eye(2);
Dxz(:,:,1,2) = Z;
Dxz(:,:,2,1) = X;
Dxz(:,:,2,2) = (-1i)*Z*X;
Dxzt(:,:,1,1) = Dxz(:,:,1,1)';
Dxzt(:,:,1,2) = Dxz(:,:,1,2)';
Dxzt(:,:,2,1) = Dxz(:,:,2,1)';
Dxzt(:,:,2,2) = Dxz(:,:,2,2)';
A00 = 1/2*(Dxz(:,:,1,1) + Dxz(:,:,1,2) + Dxz(:,:,2,1) + Dxz(:,:,2,2));
A = zeros(2,2,2,2);
for i = 1:2
    for j = 1:2
        A(:,:,i,j) = Dxz(:,:,i,j)*A00*Dxzt(:,:,i,j);
    end
end
for i = 1:2
    for j = 1:2
        M = 0;
        for k = 1:2
            for l = 1:2
                M = M + 1/2*abs(trace(A(:,:,k,l)*U*A(:,:,i,j)*Ut));
            end
        end
        if mana(1) < M
            mana(1) = M;
        end
    end
end

%---d = 3------------------------------------------------------------------
U = diag([1, exp(1i*2*pi/9), exp(1i*16*pi/9)]);
Ut = U'; % ' automatically implies conjugate transpose
inv2 = modinv(2,3);
X = [zeros(1, 2) 1; eye(2, 3)];
Z = diag(exp(1i*2*pi/3).^(0:2));
Dxz = zeros(3,3,3,3);
Dxzt = zeros(3,3,3,3);
A00 = zeros(3,3);
for x = 0:2
    for z = 0:2
        Dxz(:,:,x+1,z+1) = exp(-1i*2*pi/3*inv2*x*z)*(Z^z)*(X^x);
        Dxzt(:,:,x+1,z+1) = Dxz(:,:,x+1,z+1)';
        A00 = A00 + 1/3*Dxz(:,:,x+1,z+1);
    end
end
A = zeros(3,3,3,3);
for i = 1:3
    for j = 1:3
        A(:,:,i,j) = Dxz(:,:,i,j)*A00*Dxzt(:,:,i,j);
    end
end
for i = 1:3
    for j = 1:3
        M = 0;
        for k = 1:3
            for l = 1:3
                M = M + 1/3*abs(trace(A(:,:,k,l)*U*A(:,:,i,j)*Ut));
            end
        end
        if mana(2) < M
            mana(2) = M;
        end
    end
end

%---d = n------------------------------------------------------------------
for n = 3:np;
    d = p(n);
    U = diag(exp(1i*2*pi/d*((0:d-1).^3)));
    Ut = U';
    inv2 = modinv(2,d);
    X = [zeros(1, d-1) 1; eye(d-1, d)];
    Z = diag(exp(1i*2*pi/d).^(0:d-1));
    Dxz = zeros(d,d,d,d);
    Dxzt = zeros(d,d,d,d);
    A00 = zeros(d,d);
    for x = 0:d-1
        for z = 0:d-1
            Dxz(:,:,x+1,z+1) = exp(-1i*2*pi/d*inv2*x*z)*(Z^z)*(X^x);
            Dxzt(:,:,x+1,z+1) = Dxz(:,:,x+1,z+1)';
            A00 = A00 + 1/d*Dxz(:,:,x+1,z+1);
        end
    end
    A = zeros(d,d,d,d);
    for i = 1:d
        for j = 1:d
            A(:,:,i,j) = Dxz(:,:,i,j)*A00*Dxzt(:,:,i,j);
        end
    end
    for i = 1:d
        for j = 1:d
            M = 0;
            for k = 1:3
                for l = 1:d
                    M = M + 1/d*abs(trace(A(:,:,k,l)*U*A(:,:,i,j)*Ut));
                end
            end
            if mana(n) < M
                mana(n) = M;
            end
        end
    end
end

figure(1)
labels = strcat(repmat('(',np,1), num2str(p'), repmat(', ',np,1), num2str(mana'), repmat(')',np,1));
plot(p,mana,'--gs','LineWidth',2,'MarkerSize',10,'MarkerEdgeColor','b','MarkerFaceColor',[0.5,0.5,0.5])
text(p, mana, labels, 'VerticalAlignment','top', 'HorizontalAlignment','left');