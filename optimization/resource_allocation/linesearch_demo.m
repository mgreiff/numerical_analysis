%% x=[0;0] and d=[1;0]
disp('~~~ test_func tests using Armijjo~~~')
disp('For x=[0;0] and d=[1;0] using test_func...')
x=[0;0];
d=[1;0];
tic
[lambda, numiter] = linesearchArmijjo(@test_func , x , d)
fval = test_func(x + lambda.*d)
disp(['Was computed in ', num2str(toc) , ' seconds'])

%% x=[0;0] and d=[0;1]
disp('For x=[0;0] and d=[0;1] using test_func...')
x=[0;0];
d=[0;1];
tic
[lambda, numiter] = linesearchArmijjo(@test_func , x , d)
fval = test_func(x + lambda.*d)
disp(['Was computed in ', num2str(toc) , ' seconds'])

%% x=[0;0] and d=[1;0]
disp('~~~ test_func tests using GS~~~')
disp('For x=[0;0] and d=[1;0] using test_func...')
x=[0;0];
d=[1;0];
tic
[lambda, numiter] = linesearchGS(@test_func , x , d)
fval = test_func(x + lambda.*d)
disp(['Was computed in ', num2str(toc) , ' seconds'])

%% x=[0;0] and d=[0;1]
disp('For x=[0;0] and d=[0;1] using test_func...')
x=[0;0];
d=[0;1];
tic
[lambda, numiter] = linesearchGS(@test_func , x , d)
fval = test_func(x + lambda.*d)
disp(['Was computed in ', num2str(toc) , ' seconds'])

%% x=[0;0] and d=[1;0]
disp('~~~ test_func tests using Newton~~~')
disp('For x=[0;0] and d=[1;0] using test_func...')
x=[0;0];
d=[1;0];
tic
[lambda, numiter] = linesearchNewton(@test_func , x , d)
fval = test_func(x + lambda.*d)
disp(['Was computed in ', num2str(toc) , ' seconds'])

%% x=[0;0] and d=[0;1]
disp('For x=[0;0] and d=[0;1] using test_func...')
x=[0;0];
d=[0;1];
tic
[lambda, numiter] = linesearchNewton(@test_func , x , d)
fval = test_func(x + lambda.*d)
disp(['Was computed in ', num2str(toc) , ' seconds'])

%~~~~~~ Final linesearch implementation ~~~~~~~
%% x=[0;0] and d=[1;0]
disp('~~~ test_func tests using Newton~~~')
disp('For x=[0;0] and d=[1;0] using test_func...')
x=[0;0];
d=[1;0];
tic
[lambda, numiter] = linesearch(@test_func , x , d)
fval = test_func(x + lambda.*d)
disp(['Was computed in ', num2str(toc) , ' seconds'])

%% x=[0;0] and d=[0;1]
disp('For x=[0;0] and d=[0;1] using test_func...')
x=[0;0];
d=[0;1];
tic
[lambda, numiter] = linesearch(@test_func , x , d)
fval = test_func(x + lambda.*d)
disp(['Was computed in ', num2str(toc) , ' seconds'])
%% x=[0;0] and d=[0;1]
disp('For x=[0] and d=[1] using test_func2...')
x=[0];
d=[1];
f=@(x)(1-1e-10*x)^2;
tic
[lambda, numiter] = linesearch(f , x , d)
fval = f(x + lambda.*d)
disp(['Was computed in ', num2str(toc) , ' seconds'])