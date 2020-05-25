function f = deformedschaffer2(x)

temp1 =x(:,1:end-1).^2 - x(:,2:end).^2;
temp2 = x(:,1:end-1).^2 + x(:,2:end).^2;

fact1 = (sin(temp1)).^2 - 0.5;
% fact2 = (1 + 0.001*(temp2)).^2;
fact2 = (1 + 0.01*(temp2)).^2; % change the coefficient value
f= 0.5 + sum( fact1 ./ fact2,2);
return