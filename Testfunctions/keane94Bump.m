function f = keane94Bump (X)

[s1,d] = size(X);

f1 = sum(cos(X).^4,2);
f2 = prod(cos(X).^2,2);
f3 = sum(repmat(1:d,s1,1).* (X.^2),2).^0.5;

g1 = 0.75 - prod(X,2);
g2 = sum(X,2) - 7.5.*d;

temp_a = g1 > 0 ;
temp_b = g2 > 0;


 f_original = -abs((f1 - 2.*f2) ./ f3);
 f = f_original +  temp_a.* (g1.* f_original).^2 + temp_b.* (g2 .* f_original).^2;
 
 return


