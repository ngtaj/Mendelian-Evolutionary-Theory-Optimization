function [y] = cont(X, u, a)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% CONTINUOUS INTEGRAND FAMILY
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% INPUTS:
%
% xx = [x1, x2, ..., xd]
% u = [u1, u2, ..., ud] (optional), with default value
%     [0.5, 0.5, ..., 0.5]
% a = [a1, a2, ..., ad] (optional), with default value [5, 5, ..., 5]
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


[s1, d] = size(X);

if (nargin == 1)
    u = repmat(0.5, s1, d);
    a = repmat(5, s1, d);
elseif (nargin == 2)
    a = repmat(5, 1, d);
end
y = exp(-sum(a .* abs(X-u),2));
% y=zeros(s1,1);
% for i = 1: s1
%     
%     xx = X(i,:);
% sum = 0;
% for j = 1:d
%    xi = xx(j);
%    ai = a(j);
%    ui = u(j);
%    new = ai .* abs(xi-ui);
%    sum = sum + new;
% end
% 
% y(i,1) = exp(-sum);
% 
% end
return
