function f = ackley(X)


%               -35 < x_i < 35
n= size(X,2);
            f = 20.*(1 - exp(-0.2*sqrt((1/n).*sum(X.^2 , 2))))...
                   - exp((1/n).*sum(cos((2*pi).*X) , 2)) + exp(1);
  
return