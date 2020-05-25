function f = rosenbrock_Leon(x)
% Extended Rosenbruck's Banana-function, for N-dimensional input
%
%   ROSENBRUCK([x1, x2, .., xn]) returns the value of the Rosenbruck
%   function at the specified points. All [xi] may be vectors. The search 
%   domain is
%
%               -100 < x_i < 100
%
%   The global minimum is 
%
%               f(x1, x2, ..., xn) = f(1, 1, ..., 1) = 0


        f = sum(100.*(x(:,2:end) - x(:,1:end-1).^3).^2 + (1 - x(:,1:end-1)).^2,2); 


return