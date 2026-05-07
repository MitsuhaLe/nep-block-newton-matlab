% example 15
% T(lambda) = A_1 + f(lambda)* e_n e_n' - lambda * A_3
clear; clc

% Set the problem 
n = 100;
e = ones(n,1);
A1 = spdiags([-e 2*e -e],-1:1,n,n);
A1(n,n) = 1;
A1 = n* A1;
A3 = spdiags([e 4*e e],-1:1,n,n);
A3(n,n) = 2;
A3 = 1/(6*n) * A3;
A2 = sparse(n,n);
A2(n,n) = 1;
A = zeros(n,n,3);
A(:,:,1) = A1;
A(:,:,2) = A2;
A(:,:,3) = A3;

% Initialization
k = 5; l = 1;
S = diag(2*ones(k,1));
% 
% X = randn(n,k);
% 
% W = zeros(n,k,l);
% for i = 1:l
%     W(:,:,i) = X * (S^(i-1));
% end
%
[W, X, S, res] = initialPair(A, @f, S, l);
WW = zeros(n,k,l);
for i = 1:l
    WW(:,:,i) = W((i-1)*n + (1:n),:);
end

%%% Main Alg
iterMax = 20;
tol = 1e-10;
Res = zeros(iterMax+1,1);
stepLength = zeros(iterMax+1,1);
eigVals = zeros(iterMax+1, k);

for iter = 1:iterMax
    [RT, res] = resT(A, @f, X, S);
    fprintf('Iteration %d, residual = %.5e\n', iter-1, res);
    Res(iter) = res;
    RV = zeros(k);
    [U, RR] = schur(S);
    eigVals(iter,:) = diag(RR)';
    if  res < tol
        break
    end
    [dX, dS] = nlevp_newtonstep( A, @f, X*U, RR, W, RT*U, RV );
    dX = dX * U';
    dS = U * dS * U';
    % use armijo rule
    [j, tau] = armijoRule(A, @f, X, S, dX, dS, res);
    fprintf("step_length = %g\n", tau);
    stepLength(iter) = tau;
    X = X - dX * tau;
    S = S - dS * tau;
    % Orthogonalization progress
    V = zeros(n*l,k);
    V(1:n,:) = X;
    for i = 1:l-1
        V(i*n + (1:n), :) =  V((i-1)*n + (1:n),:) * S;
    end
    [W, R] = qr(V,"econ");
    V = zeros(n,k,l);
    for i = 1:l
        V(:,:,i) = W((i-1)*n + (1:n),:);
    end
    W = V;
    X = X/R;
    S = (R * S)/R;
end
fprintf('final iteration %d, residual = %.5e\n', iter, res);

figure(1)
semilogy(Res)
xlabel("iteration")
ylabel("Residual")

figure(2)
eigVals = sort(eigVals(1:iter,:), 2);
plot(eigVals, 'LineWidth', 1.2)
xlabel("iteration")
ylabel("Eigenvalue")
grid on




function y = f(j,S)
% M = K = kappa = 1
    n = size(S,1);
    switch j
        case 1
            y = eye(n);
        case 2
            y = eye(n) + (S - eye(n))^(-1);
        case 3
            y = -S;
    end
end