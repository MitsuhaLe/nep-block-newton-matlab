function [W, X, S, res] = initialPair(A, f, S, l)
% algorithm 2 for computing initial X
%   S is supposed to be upper triangular matrix
n = size(A,1);
m = size(A,3);
k = size(S, 1);
rng(12)
X = randn(n,k);
Y = zeros(n,k);


res = zeros(4,0);
[~, res(1)] = resT(A, f, X, S);
for p = 1:3
    % compute f(j,S)
    fS = zeros(k,k,m);
    for j = 1:m 
        fS(:,:,j) = feval(f,j,S); 
    end
    % step 3
    for i = 1:k
        s = S(i,i);
        coeffMatrix = zeros(n);
        for j = 1:m
            coeffMatrix = coeffMatrix + A(:,:,j)*feval(f, j, s); 
        end
        rhs = X(:,1);
        Y(:,i) = coeffMatrix \ rhs;
        Z = zeros(n,k-i+1);
        Z(:,1) = Y(:,i);
        for j = 1:m
            X = X - A(:,:,j) * Z * fS(i:end,i:end,j);
        end
        X = X(:, 2:end);
    end
    % step 4,5
    % Orthogonalization progress
    V = zeros(n*l,k);
    V(1:n,:) = Y;
    for i = 1:l-1
        V(i*n + (1:n), :) =  V((i-1)*n + (1:n),:) * S;
    end
    [W, R] = qr(V,"econ");
    X = Y/R;
    S = (R * S)/R;
    [~, res(p+1)] = resT(A, f, X, S);
end

end

