%% Solve auxillary problem
f = @(x) exp(x(1).*x(2).*x(3).*x(4).*x(5));

h1 = @(x) x(1).^2+x(2).^2+x(3).^2+x(4).^2+x(5).^2-10;
h2 = @(x) x(2).*x(3) - 5.*x(4).*x(5);
h3 = @(x) x(1).^3 + x(3).^3 + 1;

alpha = @(x) (h1(x).^2 + h2(x).^2 + h3(x).^2);

% ~~~ solver Options ~~~
method = 'BFGS';           % Method (DFP,BFGS,Newton,fminunc)
tol = 1e-4;                   % Tolerance
x = [-2;2;2;-1;-1];               % Initial point
mu = [0.01,0.1,1,10,100];     % Strategy
printout = 0;
solhist = [x];
m = 0.01*2.*[h1(x);h2(x);h3(x)];
functionValues = [f(x)];

% Solver
for ii =1:length(mu)
    aux = @(x) f(x) + mu(ii).*alpha(x);
    
    if strcmp(method, 'fminunc')
        x = fminunc(aux , x);
    else
        [x, ~] = nonlinearmin(aux , x , method , tol , printout);
    end
    
    solhist = [solhist,x];
    functionValues = [functionValues f(x)];
end

h0 = [h1(solhist(:,1));h2(solhist(:,1));h3(solhist(:,1))];
multipliers = mu(end)*2.*[h1(solhist(:,end));h2(solhist(:,end));h3(solhist(:,end))];