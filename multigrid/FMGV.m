function [v , state] = FMGV( v , f , currentGrid, state)

state.gridHistory = [state.gridHistory,currentGrid];

N = length(v)-2;
h = 1/(N + 1);

% ~~~ Uncomment to use the gamma sequence ~~~
% indlength = 20;
% gammasec = linspace(0.95,0.6,indlength);
% indices = 3:(indlength+2);
% if N < 2^3-1
%     gamma = 0.95;
% elseif N > 2^indices(end)-1
%     gamma = 0.6;
% else
%     for ii = 1:length(indices)
%         if N == 2^indices(ii) - 1
%             gamma = gammasec(ii);
%         end
%     end
% end
gamma = 2/3;

if state.numberOfGrids == currentGrid
    I = speye(N,N);
    subdiag = sparse(2:N,1:N-1,1,N,N);
    D = I/2-(state.beta/h^2)*(subdiag+subdiag'-4*I/2);
    Teff = kron(D,I)+kron(I,D);
    fhat = reshape(f(2:end-1,2:end-1),N^2,1);
    v(2:end-1,2:end-1) = reshape(Teff\fhat,N,N);
    return; 
end

kernel = state.beta*[0,-1,0;-1,4,-1;0,-1,0]/h.^2;
kernel(2,2) = kernel(2,2) + 1;



for ii = 1:state.eta1
    
    rf = v;

    % Compute residual
    rf(2:end-1,2:end-1) = conv2(rf,kernel,'valid') - f(2:end-1,2:end-1);

    % Pre-smoothing Jacobi
    v(2:end-1,2:end-1) = v(2:end-1,2:end-1) - gamma*(1/kernel(2,2))*rf(2:end-1,2:end-1);
end

% Compute residual
rf(2:end-1,2:end-1) = conv2(v,kernel,'valid') - f(2:end-1,2:end-1);


% Anti-alias filter
rf = lowpass(rf);

% Restrict to coarse grid
rc = FMGrestrict(rf);

for ii = 1:state.nVW
% Solve error equation
[ec, state] = FMGV(0*rc , rc , currentGrid + 1 , state);
state.gridHistory = [state.gridHistory,currentGrid];
end

% Prolong to one grid
ef = FMGprolong(ec);

% Correct and remove error
v = v - ef;

for ii = 1:state.eta2
    % Compute residual
    rf(2:end-1,2:end-1) = conv2(v,kernel,'valid') - f(2:end-1,2:end-1);

    % Post-smoothing Jacobi
    v(2:end-1,2:end-1) = v(2:end-1,2:end-1) - gamma*(1/kernel(2,2))*rf(2:end-1,2:end-1);
end
end

