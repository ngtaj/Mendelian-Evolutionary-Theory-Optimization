 function [f] = schaffer2(x)



% SCHAFFER FUNCTION N. 2
%
% INPUT:
%
% xx = [x1, x2]
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[s1,d] = size(x);


fact1 = (sin(x(:,1:end-1).^2 - x(:,2:end).^2)).^2 - 0.5;
fact2 = (1 + 0.001*(x(:,1:end-1).^2+x(:,2:end).^2)).^2;

f = 0.5 + sum(fact1./fact2,2);





end

