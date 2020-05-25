function f = LunacekBiRastrigin(X)

[s1,n] = size(X);
s=1-1/(2.*sqrt(n+20)-8.2); % Or, it can be any float value between 0.2 and 1.4 as in the original reference
% s_2=1-1/(2*sqrt(2+20)-8.2); 
d=1; % 1, 2, 3 or 4
mu1=2.5;
mu2_1=-sqrt((mu1^2-d)/s); % for 1-dimensional plotting
f= zeros(s1,1);
for i = 1: s1
x = X(i,:);

f(i,1)= min([sum((x-mu1).^2,2), d.*n + s.*sum((x-mu2_1).^2)]) +...
    10.*sum((1-cos(2.*pi.*(x-mu1))),2);

end
return