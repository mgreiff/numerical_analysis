%% Plot problem and constraints
% Objective function
f = @(x) (x(1) - 5).^2 + (x(2) - 3).^2;

x1 = linspace(-2,4,100);
x2 = linspace(-2,5,100);
[X1,X2] = meshgrid(x1,x2);
for ii =1:length(x1)
    for jj =1:length(x2)
        Z(ii,jj) = f([X1(ii,jj),X2(ii,jj)]);
    end
end
figure(1)
hold on;
contour(X1,X2,Z)

% Constraint 1
cx1 = @(x) 3 - x;
xx = linspace(-2,4,100);
plot(xx,cx1(xx),'b')

% Constraint 2
cx2 = @(x) (4 + x)./2;
xx = linspace(-2,4,100);
plot(xx,cx2(xx),'r')
colorbar;

%% Solve and plot auxillary problem
% ~~~ solver Options ~~~
method = 'DFP';           % Method
tol = 1e-6;                   % Tolerance
x = [2;-1.5];                   % Initial point
epsilon = [1,.1,.01]; % Strategy
solhist = [x];

% Solver
for ii =1:length(epsilon)
    aux = @(x) f(x) + epsilon(ii).*barrierP3(x);
    
    if strcmp(method, 'fminunc')
        x = fminunc(aux , x);
    else
        [x , ~] = nonlinearmin(aux , x , method , tol , 0);
    end
    
    solhist = [solhist,x];
end
hold on;
% Visualize results
for ii =1:(length(solhist(1,:))-1)
    plot(solhist(1,ii),solhist(2,ii),'r*')
    plot([solhist(1,ii),solhist(1,ii+1)],...
         [solhist(2,ii),solhist(2,ii+1)],'r')
end

