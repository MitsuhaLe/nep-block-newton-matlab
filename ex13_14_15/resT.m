function [resMat, resNorm] = resT(A, f, X, S)
    resMat = 0*X;
    m = size(A,3);
    for i = 1:m
        resMat = resMat + A(:,:,i) * X * feval(f, i, S);
    end
    resNorm = norm(resMat, 'fro');
end
