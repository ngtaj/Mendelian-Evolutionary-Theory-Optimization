function f = holdertable(x)

% Holder table function
%
%   HOLDERTABLE([x1, x2]) returns the value of the Holder table
%   function at the specified points. [x1] and [x2] may be vectors.
%   The search domain is
%
%               -10 < x_i < 10
%
%   The four global minima are near the edges of the interval, and have a
%   function value of 
%
%       f(x*) = -1.92085026. 

% 
        % output function value
        f = -sum(abs(sin(x(:,1: end-1)).*cos(x(:,2:end)).*exp(abs(1 - sqrt(x(:,1: end-1).^2 + x(:,2:end).^2)./pi))),2);
        
