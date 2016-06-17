function [ solx , fval ] = nonlinearmin(func , x , method , tol , printout)
    %% Performs a newton iteration with various quasinewton methods to
    %  solve the unconstrained minimization of the vecotr value function.
    % ARGS
    %     func: A vector valued function handle (@(.))
    %     start: An initial position for the newton iteration (1xn array)
    %     method: Specifies which method to use when updating the inverted
    %         hessian, can be set to 'BFGS', 'DFP' or 'Newton' (str)
    %     tol: Numerical tolerance (float)
    %     printout: Specified if the function should print data on each
    %         iteration (boolean)
    % RETURNS
    %     solx: A vector solving the problem (nx1 array)
    %     fval: The function value when minimized (float)
    outerLimit = 1e3;
    deltax = 1e-8;
    solx = x;
    convHist = nan(length(x),5);
    
    if printout
        disp('    Outer iter. | Inner iter. | alpha | func(x) | stepiter\n');
    end

    % Outer iteration
    for iteration = 1:outerLimit
        
        grad = compute_gradient(func , x , deltax);
        invH = eye(length(x));

        if any(isnan(grad))
            disp('here')
        end
            
        % Inner iteration
        for ii = 1:length(x)

            s = -invH*grad; % Approximate direction
            [alpha, stepiter] = linesearch(func , x , s);
            xNew = x + alpha.* s;
            
            convHist(:,2:end)  = convHist(:,1:end-1);
            convHist(:,1) = xNew;
            
            gradNew = compute_gradient(func , xNew , deltax);
            
             % Update inverted hessian
            if strcmp(method, 'BFGS')
                invHNew = update_BFGS(invH , xNew - x , gradNew - grad);
            elseif strcmp(method, 'DFP')
                invHNew = update_DFP(invH , xNew - x , gradNew - grad);
            elseif strcmp(method, 'Newton')
                invHNew = inv(compute_hessian(func , xNew , deltax));
            else
                error(['The method "', method, '" is an invalid method.'])
            end
            
            % Prints data
            if printout
                fprintf('%12.0f %12.0f %12.10f %12.4f %12.0f\n',iteration,ii,alpha,func(x),stepiter);
            end

            % Check solution criterea
            if norm(grad) < tol
                if printout
                    disp(['Solution found in ', num2str(iteration), ' iterations with ', method])
                end
                fval = func(x);
                solx = x;
                return;
            end
            if all(~isnan(convHist(:,end))) && all(abs(convHist(:,1) - convHist(:,end)) < tol)
                if printout
                    disp(['Solution found in ', num2str(iteration), ' iterations with ', method])
                end
                fval = func(x);
                solx = x;
                return;
            end
            
            % Updates varialbles
            x = xNew;
            grad = gradNew;
            invH = invHNew;
        end
        
        solx = x;
        fval = func(solx);
    end
    if iteration == outerLimit
        disp(['Waring. Reached iteration limit without converging to a',...
             'solution with the specified tolerance.'])
        return;
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
function hess = compute_hessian(func , x , deltax)
    %% Computes the hessian of a vectorvalued function using a second order
    %  central differences scheme.
    % ARGS
    %     func: A vector valued function handle (@(.))
    %     x: A positional vector at which the gradient is approximated (nx1 array)
    %     deltax: The step size used in the central differences scheme (float)
    % RETURNS
    %     hess: The function hessian (nxn array)

    dim = length(x);
    hess = zeros(dim, dim);
    for ii = 1:dim
        for jj = 1:dim
            if jj == ii
                %Diagonal elements - pure derivatives
                hess(ii, jj) = ((func(shift(x, ii, deltax))...
                               - 2 .* func(x)...
                               + func(shift(x, ii, -deltax)))...
                               / deltax.^2);
            else
                %Other elements - mixed derivatives
                posShifted = shift(x, ii, deltax);
                negShifted = shift(x, ii, -deltax);
                hess(ii, jj) = ((func(shift(posShifted, jj, deltax))...
                               - func(shift(posShifted, jj, -deltax))...
                               - func(shift(negShifted, jj, deltax))...
                               + func(shift(negShifted, jj, -deltax)))...
                               / (4 .* deltax .^ 2));
            end
        hess = (0.5 .* (hess + hess')); % To guarantee symmetry so that the 
                                        % cholesky decomposition can be used
        end
    end
end
function hess = update_DFP(H , dX , dG)
    %% A DFP update of the inverted hessian
    % ARGS
    %     H: The old inverted hessian (nxn array)
    %     dX: The change in position (nx1 array)
    %     dG: The change in gradient (nx1 array)
    % RETURNS
    %     hess: The updated inverted hessian (nxn array)
    dxdg = (dX'*dG);
    if dxdg == 0
        dxdg = 1e-5;
    end
    dgHdg = (dG'*H*dG);
    if dgHdg == 0
        dgHdg = 1e-5;
    end
    hess = H + ((dX*dX')/dxdg) - ((H*dG*dG'*H)/dgHdg);
end
function hess = update_BFGS(H , dX , dG)
    %% A BFGS update of the inverted hessian
    % ARGS
    %     H: The old inverted hessian (nxn array)
    %     dX: The change in position (nx1 array)
    %     dG: The change in gradient (nx1 array)
    % RETURNS
    %     hess: The updated inverted hessian (nxn array)
    dxdg = (dX'*dG);
    if dxdg == 0
        dxdg = 1e-5;
    end
    hess = H + (1 + (dG'*H*dG)/(dxdg))*((dX*dX')/(dxdg))-...
           ((dX*dG'*H+H*dG*dX')/(dxdg));
end
