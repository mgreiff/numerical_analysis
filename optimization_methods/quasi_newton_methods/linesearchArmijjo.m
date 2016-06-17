function [lambda, numiter] = linesearchArmijjo(func , x , d)
    %Reference to algorithm formulation:
    % http://tosca.cs.technion.ac.il/book/samples/optimization_5_3.pdf
    
    lambda = 1e60; % Initial stepsize
    sigma = 0.3;
    beta = 0.83;
    imax = 3000;
    
    g = compute_gradient(func, x, 1e-2);
    g(isnan(g))=1;

    for numiter  = 1:imax
        if func(x + lambda.*d) > func(x) + sigma .* lambda .* g' * d
            lambda = lambda * beta;
        else
            return
        end
    end
end
function xshift = shift(x, index, deltax)
    %% Shifts an array element at a specified index by deltax
    % ARGS
    %     func: A vector valued function handle (@(.))
    %     x: A positional vector at which the gradient is approximated (nx1 array)
    %     deltax: The step size used in the central differences scheme (float)
    % RETURNS
    %     xshift: The shifted array (array)
    xshift = x;
    xshift(index) = xshift(index) + deltax;
end
function grad = compute_gradient(func , x , deltax)
    %% Computes the gradient of a vectorvalued function using a second order
    %  central differences scheme
    % ARGS
    %     func: A vector valued function handle (@(.))
    %     x: A positional vector at which the gradient is approximated (nx1 array)
    %     deltax: The step size used in the central differences scheme (float)
    % RETURNS
    %     grad: The function gradient (array)

    dim = length(x);
    grad = zeros(dim,1);
    for ii = 1:dim
        grad(ii) = (func(shift(x, ii, deltax)) - func(shift(x, ii, -deltax))) / (2 * deltax);
    end
end