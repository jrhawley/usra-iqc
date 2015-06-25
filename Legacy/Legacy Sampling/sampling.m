%---Global Declarations
global n % dimension
global w % primitive root of unity
global I % identity
global X % Weyl-Heisenberg X
global Z % Weyl-Heisenberg Z
global A % phase-point operators
global E % measurement projectors

%---Assignment
n = 3;
w = exp(1i*2*pi/n);
I = eye(n);
X = [zeros(1, n-1) 1; eye(n-1, n)];
Z = diag(w.^(0:n-1)');
A = ppo(n);
E = zeros(n,n,n); %basis measurement

[V,D] = eig(X);

for j = 1:n
    E(:,:,j) = V(:,j)*V(:,j)';
end

%---Start
% randomly generate the state
%[V, D] = eig(Z);
%c = zeros(n,1);
%for j = 1:n
%    c(j) = randi([1 n], 1);
%end
%psi = V*c;
psi = rand(n,1);
psi = psi/norm(psi);
rho = psi*psi';
U = eye(n); % orthogonal matrix U

% generate state and measurement wigner functions
[W_rho, M_rho] = dwigner(rho);
[W_E, M_E] = ewigner(E)
[W_U, M_U] = uwigner(U);
W_rho_est = zeros(n);
rho_est = zeros(n);

% correct signs for small values (floating point errors)
M_rho = correct(M_rho);
for j = 1:n
    for k = 1:n
        W_rho(j,k) = correct(W_rho(j,k));
        M_E(j,k) = correct(M_E(j,k));
        for l = 1:n
            W_E(j,k,l) = correct(W_E(j,k,l));
            for m = 1:n
                W_U(j,k,l,m) = correct(W_U(j,k,l,m));
            end
        end
    end
end

%---Sampling Protocol
T = 10; %repetitions
t = 1:T;
p_alpha = hmatcat(abs(W_rho)/M_rho); %probability of selecting a particular box
r = zeros(T,2);

% Step 2
Y = randsample(1:n^2, T, true, p_alpha); %sample from the weighted distribution
alpha = indlm1(n,Y) %convert from linear indexing to components
    
% Step 3
p_beta_alpha = zeros(T,n^2);
for j = 1:n^2
    for k = 1:T
        v = indlm1(n,j);
        p_beta_alpha(k,j) = abs(W_U(v(1), v(2), alpha(k,1), alpha(k,2)))/M_U(alpha(k,1), alpha(k,2));
    end
end
for j = 1:T
    Y(j) = randsample(1:n^2, 1, true, p_beta_alpha(j,:)); %sample from the weighted distribution
end
beta = indlm1(n,Y) %convert from linear indexing to components

% Step 4
p_measure = zeros(T,n);
for j = 1:T
    k = 1:n;
    p_measure(j,k) = abs(W_E(beta(j,1), beta(j,2), k))
    r(j,2) = randsample(1:n, 1, true, p_measure(j,:));
end

% Step 5
beta = v;
for j = 1:T
    r(j,1) = sign(W_rho(alpha(j,1), alpha(j,2))*W_U(beta(j,1), beta(j,1), alpha(j,1), alpha(j,2))*W_E(beta(j,1), beta(j,2), r(j,2)))*M_rho*M_U(alpha(j,1), alpha(j,2));
end

psi_est = zeros(n,1);
for j = 1:n
  psi_est(j) = mean(r(:,1).*kdelta(j,r(:,2)));
end
U*psi
psi_est = sqrt(abs(psi_est))
% Count outcomes given by sampling (ie messages sent from Alice to Bob)
%for j = 1:T
%    W_rho_est(alpha(j,1), alpha(j,2)) = W_rho_est(alpha(j,1), alpha(j,2)) + alpha(j,3);
%end
%W_rho_est = W_rho_est/sum(sum(W_rho_est)) %normalize matrix

% Bob's approximate recreation of rho from Wigner matrix estimate
%for j = 1:n 
%    for k = 1:n
%        rho_est = rho_est + W_rho_est(j,k)*A(:,:,j,k);
%    end
%end

% comparison of actual and approximate density matrices
%norm(rho - rho_est)
        
%Y = randsample(1:n, T, true, p_outcome(:,t)); %sample from weighted distribution
%    r = r*M_rho*sign(W_E(v(1), v(2), Y)); %update r
    
    % Step 4
%    R(v(1), v(2), Y, t) = r; % R(alpha, o, t) = r

%P = zeros(n,1);
%for o = 1:n %over all outcomes
%    for j = 1:n %over indices
%        for k = 1:n
%            for t = 1:T
%                P(o) = P(o) + R(j,k,o,t)*p_alpha(indml(n,[j;k]))*p_outcome(o,t);
%            end
%        end
%    end
%end
%psi_est = sqrt(P/sum(P));
%[psi psi_est]