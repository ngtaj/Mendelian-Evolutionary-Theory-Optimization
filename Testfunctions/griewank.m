function f = griewank(X)


% Griewank funcion 
%             -100 < x_i < 100
%
%   The global minimum is 
%
%               f(x) = f(0, 0) = 0.    
[m, n] = size(X);


        f =  sum(X.^2 , 2)./4000  - prod (cos (X./sqrt(repmat(1:n,m,1))) , 2) + 1;
        
        
        
 return