function scores = shubert4fcn(x)
% x= [-10,10]
    n = size(x, 2);
    
    scores = 0;
    for i = 1:n
        for j = 1:5
            scores = scores + j * cos(((j + 1) * x(:, i)) + j);
        end
    end
end