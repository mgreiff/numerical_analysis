function [lambda, numiter] = linesearch(func , x , d)
    lambda = 1;
    sigma = 0.23;
    beta = 0.78;
    imax = 600;
    deltax = 1e-12;
    
    g = compute_gradient(func, x, deltax);
    g(isnan(g))=1;

    fx = func(x);
    
    if func(x + lambda.*d) > (fx + sigma .* lambda .* g' * d)
    	searchdir = 1;
    else
        searchdir = -1;
    end
    for numiter  = 1:imax
        switch searchdir
            case 1
                if func(x + lambda.*d) > (fx + sigma .* lambda .* g' * d)
                    lambda = lambda * beta;
                else
                    return;
                end
            case -1
                if func(x + lambda.*d) <= (fx + sigma .* lambda .* g' * d)
                    lambda = lambda * (1/beta);
                else
                    return;
                end
        end       
    end
    disp('Warning! Linesearch ran the full number of iterations without finding a minima.')
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