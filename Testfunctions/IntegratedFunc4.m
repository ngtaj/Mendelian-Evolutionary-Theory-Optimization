function f=IntegratedFunc4(x)
% -pi < x < pi
 [m, n]= size(x);
 
f = wavy(x) + dropwave(x);

return