function scores = gramacyleefcn(x)
    [m, n] = size(x); % x =[-0.5, 2.5]
%     assert(n == 1, 'Gramacy & Lee function is only defined on a 1-D space.')
scores = zeros(m,1);
for i = 1:n
    scores = scores + (sin(10 .* pi .* x(:,i)) ./ (2 * x(:,i)) ) + ((x(i) - 1) .^ 4);
end
end

