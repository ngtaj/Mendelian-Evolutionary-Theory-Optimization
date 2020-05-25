function f = h11(x)

% x = 3 to 13

f = sum(sin(x(:,1:end-1) + x(:,2:end)) + sin((2.*x(:,1:end-1).* x(:,2:end)) ./ 3),2);

return