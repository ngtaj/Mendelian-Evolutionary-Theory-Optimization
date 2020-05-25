function f = wavy(x)

% -pi < x < pi
k =10;

f =1 - (1/size(x,2)).*(sum(cos(k.*x).*exp(-(x.^2)./2),2));
return