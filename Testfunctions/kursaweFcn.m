function f = kursaweFcn(x)

% Kursawe Pareto front test problem

%  Copyright (c) 2010, The MathWorks, Inc.
%  All rights reserved.
% -5 < x < 5

% function 1
f1 = sum (-20.*exp(-0.2.*sqrt(x(:,1:end-1).^2 + x(:,2:end).^2)),2);

% function 2

f2 = sum(abs(x).^0.8 + 5.*sin(x.^3),2);

 f = f1+f2;
return