function X = normalizeRows(X, tol)
    normX = sqrt(sum(X.^2,2));
    normX(normX < tol) = 1;
    X = X./normX;
end