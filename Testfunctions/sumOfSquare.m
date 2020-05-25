function f = sumOfSquare(X)
[s1,s2] = size(X);
i = repmat(1:s2, s1,1);

f =sum(X.^i, 2);

return