function f = GenSchwefelNo226(X)



f = 418.9829 * size(X,2)- sum(X.*sin(sqrt(abs(X))),2);

return