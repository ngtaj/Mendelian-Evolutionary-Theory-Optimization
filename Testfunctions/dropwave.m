function f = dropwave ( x )

% -5.12 <x < 5.12
%*****************************************************************************80
%
%% P11_F evaluates the objective function for problem 11.


    rsq = x(:, 1:end-1).^2 + x(:, 2:end).^2 ;

    f = -sum( ( 1 + cos ( 12 .* sqrt ( rsq ) ) ) ./ ( 0.5 .* rsq + 2.0 ),2);
