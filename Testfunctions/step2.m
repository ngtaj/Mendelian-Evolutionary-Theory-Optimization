function f = step2(X)


% INPUT:
% xx = [x1, x2, ..., xd]

  f = sum(floor((abs(X + 0.5))) .^2 ,2);
 
%  f = sum( floor(X.^2), 2); % Step function 3
return