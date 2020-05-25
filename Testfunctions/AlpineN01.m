function f = AlpineN01(x)
% -10 < X < 10

% f = sum(abs(x.*sin(x)+0.1.*x),2);
 f =  10*size(x,2) + sum((x.^2 - 25.*cos((2*pi).*x)),2);
rsq = x(:, 1:end-1).^2 + x(:, 2:end).^2 ;

  f = -50*sum( ( 1 + cos ( 15 .* sqrt ( rsq ) ) ) ./ ( 0.5 .* rsq + 2.0 ),2) +f;
end

