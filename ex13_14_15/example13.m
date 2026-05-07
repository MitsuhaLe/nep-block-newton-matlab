% example 13
% T(lambda) = lambda - A_0 - A_1 e^(-lambda*tau)
clear; clc

A(:,:,1) = [-5 1; 2 -6];
A(:,:,2) = [1 0;0 1];
A(:,:,3) = [-2 1; 4 -1];

tau = 1;

%%% compute initial pair
A0 = [-5 1; 2 -6];
A1 = [-2 1; 4 -1];
[~, eig] = polyeig(-(A0+A1), eye(2)+A1, -A1/2, A1/6, -A1/24);
figure(1)
plot(eig,'o')
hold on

% Filter out eigenvalues with real parts less than 1 
% and find the top 5 with the largest real parts from them
idx_less_than_1 = find(real(eig) < 1);
[~, relative_top_idx] = maxk(real(eig(idx_less_than_1)), 5);
idx = idx_less_than_1(relative_top_idx);
eig = maxk(eig(real(eig) < 1), 5, 'ComparisonMethod', 'real');
plot(eig)
hold off

%%% Main Alg
n = 2; k = 5; l = 3;
S = diag(eig);
[W, X, S, res] = initialPair(A, @f, S, l);
WW = zeros(n,k,l);
for i = 1:l
    WW(:,:,i) = W((i-1)*n + (1:n),:);
end

iterMax = 20;
tol = 1e-10;
Res = zeros(iterMax+1,1);
stepLength = zeros(iterMax+1,1);

for iter = 1:iterMax
    [RT, res] = resT(A, @f, X, S);
    fprintf('Iteration %d, residual = %.5e\n', iter-1, res);
%     if  res < tol
%         break
%     end
    Res(iter) = res;
    RV = zeros(k);
    [U, RR] = schur(S);
    [dX, dS] = nlevp_newtonstep( A, @f, X*U, RR, WW, RT*U, RV );
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

figure(2)
semilogy(0:iter-1, Res(1:iter))
xlabel("iteration")
ylabel("Residual")






function F = f(j,M)
    switch j
        case 1
            F = - eye(size(M,1));
        case 2
            F = M;
        case 3
            F = -expm(-M);
    end
end




% % 定义矩阵 A0
% A0 = [ -0.8498,    0.1479,    44.37;
%         0.003756, -0.2805,  -229.2;
%        -0.1754,    0.02296,  -0.3608 ];
% 
% % 定义矩阵 A1
% A1 = [ 0.28,  0,     0;
%        0,    -0.28,  0;
%        0,     0,     0 ];
% eig = polyeig(-(A0+A1), eye(3)+A1, -A1/2);
% plot(eig,'o')

