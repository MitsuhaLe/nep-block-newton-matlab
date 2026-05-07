function [j, tau] = armijoRule(A, f, X, S, dX, dS, ...
    lastErr, armijoFactor, armijoMax)
    if nargin < 8 || isempty(armijoFactor)
        armijoFactor = 0.5;
    end
    if nargin < 9 || isempty(armijoMax)
        armijoMax = 3;
    end
    j = 0; tau = 1;
    [~, res] = resT(A, f, X - dX, S - dS);
    if res <(1-10e-4 * armijoFactor) * lastErr
        return
    end
    while j < armijoMax
        j = j+1; tau = tau * armijoFactor;
        XX = X - dX * tau;
        SS = S - dS * tau;
        [~, res] = resT(A, f, XX, SS);
        if res <(1-10e-4 * tau) * lastErr
            return
        end
    end
end