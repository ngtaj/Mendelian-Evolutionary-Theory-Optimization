function f = IntegratedFunc10(x)
ps= size(x,1);
f = step2(x)+ 30*InvertedCosineWave(x);
return%
