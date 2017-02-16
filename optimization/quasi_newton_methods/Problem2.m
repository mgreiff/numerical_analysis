%% Plot problem and constraints
% Objective function
f = @(x) exp(x(1)) + x(1).^2 + x(1).*x(2);

x1 = linspace(-15,5,100);
x2 = linspace(-2,10,100);
[X1,X2] = meshgrid(x1,x2);
for ii =1:length(x1)
    for jj =1:length(x2)
        Z(ii,jj) = f([X1(ii,jj),X2(ii,jj)]);
    end
end
figure(1)
hold on;
contour(X1,X2,Z,20)

% Constraint 1
cx1 = @(x) 1-0.5.*x;
xx = linspace(-15,5,100);
plot(xx,cx1(xx),'b')

%% Solve and plot auxillary problem
alpha = @(x) (0.5.*x(1) + x(2) - 1).^2;

% ~~~ solver Options ~~~
method = 'BFGS';              % Method (DFP,BFGS,Newton,fminunc)
tol = 1e-4;                   % Tolerance
x = [4;4];                    % Initial point
mu = [4,40,400];              % Strategy
solhist = [x];

% Solver
tic
for ii =1:length(mu)
    aux = @(x) f(x) + mu(ii).*alpha(x);
    
    if strcmp(method, 'fminunc')
        x = fminunc(aux , x);
    else
        [x, ~] = nonlinearmin(aux , x , method , tol , 0);
    end
    
    solhist = [solhist,x];
end
toc

% Visualize results
for ii =1:(length(solhist(1,:))-1)
    plot([solhist(1,ii),solhist(1,ii+1)],...
         [solhist(2,ii),solhist(2,ii+1)],'g')
end
plot(solhist(1,1),solhist(2,1),'g+')
plot(solhist(1,end),solhist(2,end),'g+')

fopt = f(solhist(:,end))
lagrangeMultiplier = 8.*(0.5.*solhist(1,end) + solhist(2,end) - 1)