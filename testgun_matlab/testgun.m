% Reproduce the Nep-Pack (Julia) example: 
% nep = nep_gallery("nlevp_native_gun")
%
% A short introduction for NEP "gun" in nlevp 
%    "https://github.com/ftisseur/nlevp"
% gun: NEP from model of a radio-frequency gun cavity.
%   [COEFFS,FUN,F,XCOEFFS] = nlevp('gun') returns 9956-by-9956 matrices
%   K, M, W1 and W2 for the nonlinear eigenvalue problem defined by
%   T(lam)*x = [K-lam*M+i*lam^(1/2)*W1+i*(lam-108.8774^2)^(1/2)*W2]*x = 0.
%   The matrices are returned in a cell array: COEFFS = {K,M,W1,W2}.
%   FUN is a function handle to evaluate the functions 1,-lam,i*lam^(1/2),
%   and i*(lam-108.8774^2)^(1/2) and their derivatives.
%   F is the function handle T(lam).
%   XCOEFFS is the cell {1, 1, L1, L2; K, M, R1', R2'} to exploit the low
%   rank of W1 = L1*R1' and W2 = L2*R2'.
%   This problem has the properties nep, sparse, banded (843 bands),
%   low-rank.


% Load coefficient matrices K, M, W1, W2 
[COEFFS,~,~,~] = nlevp('gun');
m = 4;
n = size(COEFFS{1},1);
A = COEFFS;

% Initialization for S and X which come from example in nep-pack
S=150^2*[1 0; 0 1]; X=[[1 0; 0 1]; zeros(n-2,2)];
iterMax = 15;
tol = 1e-12;
Res = zeros(iterMax+1,1);
stepLength = zeros(iterMax+1,1);

% Constuct orth matrix W
k=2; l=2;
W = zeros(n,k,l);
W(:,:,1) = X;
for i = 2:l
    W(:,:,i) = W(:,:,i-1) * S;
end

% Main loop
for iter = 1:iterMax
    [RT, res] = resT(A, @f, X, S);
    fprintf('Iteration %d, residual = %.5e\n', iter-1, res);
    if  res < tol
        break
    end
    Res(iter) = res;
    RV = zeros(k);
    [U, RR] = schur(S);
    [dX, dS] = nlevp_newtonstep( A, @f, X*U, RR, W, RT*U, RV );
    % restore dx and ds
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

% define f_j(lam) in nlevp('gun') 
function y = f(j, lam)
    switch j
        case 1
            y = eye(size(lam,1));
        case 2
            y = -lam;
        case 3
            y = 1i*sqrtm(lam);
        case 4
            y = 1i*sqrtm(lam-108.8774^2 * eye(size(lam,1)));
    end
end



