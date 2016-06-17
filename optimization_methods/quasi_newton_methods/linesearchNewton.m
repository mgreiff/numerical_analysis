function [lambda, iteration] = linesearchNewton(f , x , d)

    lambdaold = 0;
    F = @(lambda) f(x + lambda.*d);
    
    delta = 1e-58;
    limit = 100;
    tol = 1e-58;
    
    for iteration = 1:limit
        dF = (F(lambdaold + delta) - F(lambdaold - delta))/(2*delta);
        ddF = (F(lambdaold + delta) -2*F(lambdaold)+ F(lambdaold - delta))/(delta.^2);
        lambda = lambdaold - dF/ddF;
        if abs(dF) < tol || abs(lambda - lambdaold) < tol
            return
        end
        lambdaold = lambda;
    end
end