function f = bird(x)
% Bird function
%
%   BIRD([x1, x2]) returns the value of the Bird function 
%   at the specified points. [x1] and [x2] may be vectors. 
%   The search domain is
%
%               -10 < x_i < 10
%
       
        % output function value
f = sum(sin(x(:,1:end-1)).*exp((1-cos(x(:,2:end))).^2) + cos(x(:,2:end)).*exp((1-sin(x(:,1:end-1))).^2) + (x(:,1:end-1)-x(:,2:end)).^2,2);
% f = prod(sin(x(:,1:end-1)).*exp((1-cos(x(:,2:end))).^2) + cos(x(:,2:end)).*exp((1-sin(x(:,1:end-1))).^2) + (x(:,1:end-1)-x(:,2:end)).^2,2);
return
