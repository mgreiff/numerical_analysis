function out = FMGprolong( v )
out = zeros(2*size(v)-1);
out(1:2:end,1:2:end) = v;
kernel = (1/16)*[1,2,1;2,4,2;1,2,1];
out = 2^2*conv2(out,kernel,'same');
end

