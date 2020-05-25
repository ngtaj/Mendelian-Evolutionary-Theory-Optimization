function f = Noisy_SinusoidalFunction(X)

% -100 to 100
A=2.5;
B=5;
Z=30;
s1 = size(X,1);
f = -( A.*prod(sind(X-Z),2)-prod(sind(B*(X-Z)),2));
f = awgn(f, 10);
return