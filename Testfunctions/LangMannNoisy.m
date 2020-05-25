function f = LangMannNoisy(x)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% LANGERMANN FUNCTION
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[s1,d]=size(x);
  m = 5;
 c = [1; 2; 5; 2; 3];
 ai = [3;5;2;1;7];
 ai_1=[5;2;1;4;9];

f = zeros(s1,1);
 for i = 1:s1
     temp = (repmat(x(i,1:end-1),m,1) - repmat(ai,1,d-1)).^2 + (repmat(x(i,2:end),m,1) - repmat(ai_1,1,d-1)).^2;
   f(i,:)  =- sum(sum((c.*cos(pi.*(temp))./(exp(temp)./pi)),1),2);
   f(i,:) = awgn(f(i,:), 10);
 end 
