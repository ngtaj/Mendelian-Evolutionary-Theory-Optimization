function [y] = michal_Ext(X)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% MICHALEWICZ FUNCTION
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% INPUTS:
%
% xx = [x1, x2]
% m = constant (optional), with default value 10
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%xi ? [0, ?]


    m = 10;

[s1, d] = size(X);

y = ones(s1,1);
for i = 1: s1
xx = X(i,:);
sum = 0;

for ii = 1:d
	xi = xx(ii);
	new = sin(xi) .* (sin(ii.*xi.^2/pi)).^(2.*m);
	sum  = sum + new;
end

y(i,1) = -sum;

end
return
