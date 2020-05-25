function f = stochasticf(x)

% -5 < x< 5%
[s1,s2]=size(x);
f =sum(randn(s1,s2).*abs(x-(1./repmat(1:s2,s1,1))),2);
return