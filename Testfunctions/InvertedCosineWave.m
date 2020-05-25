function f = InvertedCosineWave(X)



f = - sum(exp(-(X(:,1:end-1).^2 + X(:,2:end).^2 + 0.5.*X(:,1:end-1).*X(:,2:end))...
    ./8).*cos(4.*sqrt(X(:,1: end-1).^2 + X(:,2:end).^2 + 0.5.*X(:,1:end-1).*X(:,2:end))),2);

return