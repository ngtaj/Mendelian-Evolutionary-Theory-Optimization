function f = LorF2(X)


[s1,s2] = size(X);
L1=repmat(5.1,s1,s2);
L2=repmat(0.5,s1,s2);
L3=repmat(4*log(2),s1,s2);
L4=repmat(0.066832364099628,s1,s2); %0.0667;
L5=repmat(0.64,s1,s2);
k=6;

f = - prod(sin(L1.*pi.*X + L2).^k.*exp(-L3.*((X-L4)./L5).^2), 2);

return