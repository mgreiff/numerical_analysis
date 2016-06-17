function out = lowpass( v )
ker = [1,2,1;2,4,2;1,2,1]/16;
out = v;
out(2:end-1,2:end-1) = conv2(v,ker,'valid');
end

