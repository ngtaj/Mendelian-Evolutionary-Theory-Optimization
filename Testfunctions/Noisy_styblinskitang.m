function f = Noisy_styblinskitang(X)

% Styblinski-Tang function
%
%   STYBLINSKYTANG([x1, x2, ..., xn]) returns the value of the 
%   Styblinski-Tang at the specified points. All [xi] may be vectors. 
%   The search domain is 
%
%               -5 < x_i < 5
%
%   The global minimum is 
%
%               f(x1, x2, ..., xn) = 
%               f(-2.903534, -2.903534, ..., -2.903534) = -39.16599 * n


% s1 = size(X,1);
        f = sum(X.^4 - 16.*X.^2 + 5.*X, 2)./2;
        
f = awgn(f, 5);

return