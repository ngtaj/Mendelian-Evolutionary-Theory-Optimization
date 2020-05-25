function f = crossintray(x)
% Cross-in-tray function
%
%   CROSSINTRAY([x1, x2]) returns the value of thecross-in-tray 
%   function at the specified points. [x1] and [x2] may be vectors. 
%   The search domain is
%
%               -10 < x_i < 10
%
%   The global minimum is found on the planes x = 0 and y = 0, with
%
%                   f(x1, x2) = f(�1.34951, �1.34951) = -2.062612.
      
        % output function value
f = - sum(0001.*(abs(sin(x(:,1: end-1)).*sin(x(:,2:end)).*exp(abs(10 - sqrt(x(:,1:end-1).^2 + x(:,2:end).^2)./pi))) + 1).^(0.1),2);

return

