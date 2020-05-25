function f = Noisy_deceptive(X)



[ s1, m]=size(X);
  f = zeros ( s1, 1 );

  beta = 2.0;
%
%  I'm just choosing ALPHA in [0,1] arbitrarily.
alpha = zeros(1,m);
  for i = 1 : m
    alpha(i) = i ./ ( m + 1 );
  end

  for i = 1 : s1
x = X(i,:);
    
f1 = 0.0;

    for j = 1 : m

      if sum( x(j) <= 0.0 )
        g = x(j);
      elseif ( x(j) <= 0.8 * alpha(j) )
        g = 0.8 - x(j) / alpha(j);
      elseif ( x(j) <= alpha(j) )
        g = 5.0 * x(j) / alpha(j) - 4.0;
      elseif ( x(j) <= ( 1.0 + 4.0 * alpha(j) ) / 5.0 )
        g = 1.0 + 5.0 * ( x(j) - alpha(j) ) / ( alpha(j) - 1.0 );
      elseif ( x(j) <= 1.0 )
        g = 0.8 + ( x(j) - 1.0 ) / ( 1.0 - alpha(j) );
      else
        g = x(j) - 1.0;
      end

      f1 = f1 + g;

    end

    f2 = f1 ./ m;
    f(i,1) = - ( f2.^beta );

  end
   f = awgn(f, 5);
return;
