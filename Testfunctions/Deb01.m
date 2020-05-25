function f = Deb01(X)

d = size(X,2);

f = - 1/d .* sum((sin(5.*pi.*X)).^6 , 2);

end