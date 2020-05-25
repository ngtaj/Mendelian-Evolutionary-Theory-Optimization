function [f]= DeflectedCorrugatedSpring(X)
% x =[0,10]
[s1,s2] = size(X);
alpha=5;
k=5;
    
f = 0.1* sum( (X - alpha).^2 - cos(k * sqrt(sum( (X - alpha).^2 ,2)))  ,2);
return

% DeVilliersGlasser02
% CrossLegTable
% XinSheYang03
% SineEnvelope
% Whitley
% Zimmerman
% Griewank
% Trefethen
% Bukin06
% Cola
% CrownedCross
% RosenbrockModified
% Stochastic