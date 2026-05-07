function [resMat, resNorm] = resT(A, f, X, S)
% resMat:  residual matrix of equation (23)
% resNorm: 2-norm residual
    resMat = 0*X;
    m = size(A,2);
    for i = 1:m
        resMat = resMat + A{i} * X * feval(f, i, S);
    end
%     resNorm = norm(resMat, 'fro');
    resNorm = norm(resMat);
end
