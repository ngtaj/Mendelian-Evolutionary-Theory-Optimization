function f = PrG3f(x)



% Search space 0< x <1
% The global minima: x* =  (1/n0.5, �, 1/n0.5), f(x*) = 1.
n = size(x,2);
f_original = - sqrt(n).^n.*prod(x,2);

h1 = sum(x.^2,2) - 1 ;
temp_a = h1 < -1e-4 | h1 > 1e-4  ;
% temp_b = h1 >= -1e-4 ;

f = f_original + 100.* temp_a.* (h1.* f_original).^2;
return

