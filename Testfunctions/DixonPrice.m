function f = DixonPrice(X)



[s1,d] = size(X);

SUM = sum(repmat(2:d,s1,1).*(2.*X(:,1:end-1).^2 - X(:,2:end)).^2, 2);

f=(X(:,1)-1).^2 + SUM;

return