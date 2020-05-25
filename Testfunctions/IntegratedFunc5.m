function f = IntegratedFunc5(x)

% f1 = schaffer2
% f2 = styblinskitang

f1 = schaffer6(x);

temp1 =x(:,1:end-1).^2 - x(:,2:end).^2;
temp2 = x(:,1:end-1).^2 + x(:,2:end).^2;

fact1 = (sin(temp1).^2 - 0.5);
  fact2 = (1 + 0.01*(temp2)).^2;
f2= 0.5 + sum( fact1./fact2,2);
 f = f1+f2;
return