function scores = periodicfcn(x)
% -2 < x < 2
    sin2x = sin(x) .^ 2;
    sumx2 = sum(x .^2, 2);
    scores = 1 + sum(sin2x, 2) -0.1 * exp(-sumx2);
    
end