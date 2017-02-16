function [ step , numiter ] = linesearchGS( func , x , d )
    
    F = @(y) func(x + y.*d);
    inter = [0,1e59];
    
    itermax = 560;
    
    alpha = (sqrt(5) - 1)/2;
    for numiter = 1:itermax
        lambda = inter(:,1) + (1 - alpha)*(inter(:,2) - inter(:,1));
        mu = inter(:,1) + alpha*(inter(:,2) - inter(:,1));
        if F(lambda) > F(mu)
            inter(:,1) = lambda;
        else
            inter(:,2) = mu;
        end
    end
    step = (inter(2) + inter(1))/2;
end