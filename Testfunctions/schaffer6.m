function [f] = schaffer6(X)



% [s1,d] = size(X);
temp = X(:,1:end-1).^2 + X(:,2:end).^2;

f = 0.5 +sum( ((sin(sqrt (temp))).^2 -0.5)./(1 + 0.001.*temp).^2,2);

