function func = test_func(x)

% To test line search for f(t)=test_func(x+t*d).
%
% Test 1: use x=[0;0], d=[1;0] 
% Test 2: use x=[0;0], d=[0;1]
%
% For both tests we have f(0)=0 and, min=-1, 
% so a reasonable functional value should be 
% between [-1,0[ (the closer to -1 the better).

func=(1e58*x(1)-1)^100+(1e-58*x(2)-1)^2-2;