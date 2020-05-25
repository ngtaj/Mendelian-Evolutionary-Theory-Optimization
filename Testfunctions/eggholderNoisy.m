function f = eggholderNoisy(x)
% generalized egg holder function
%
%   EGGHOLDER([x1, x2, ..., xn]) returns the value of the generalized 
%   eggholder function at the specified points. All [xi] may be vectors. 
%   The search domain is
%
%               -512 < x_i < 512
%
%   global minimum (for 2 variables) is at
%
%       f(x1, x2) = f(512, 404.2319) = 959.64


% Please report bugs and inquiries to: 
%
% Name       : Rody P.S. Oldenhuis
% E-mail     : oldenhuis@gmail.com    (personal)
%              oldenhuis@luxspace.lu  (professional)
% Affiliation: LuxSpace s�rl
% Licence    : BSD



% [s1,s2] = size(x);

% Id =1;
% for i =1:s1
% while Id < s2
    
f= sum(-(x(:,2:end)+47).*sin(sqrt(abs(0.5.*x(:,1:end-1)+(x(:,2:end)+47))))...
 -x(:,1:end-1).*sin(sqrt(abs((x(:,1:end-1)-(x(:,2:end)-47))))),2);

% end
% end
f= awgn(f, 5);
return       

            