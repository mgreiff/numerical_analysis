N = 2^8-1;                         % Number of internal grid points                         % Time step
h = 1/(N + 1);                     % Spatial step
[X,Y]=meshgrid(linspace(0,1,N+2)); % Mesh of the square [0,1]x[0,1]

u0 = zeros(N+2,N+2);
f = sin(2 * pi * X).*sin(pi*Y);
dt = 0.01;
beta = dt^2/4;

% Enforces boundary conditions
u0(:,1) = 0; u0(:,end) = 0; u0(1,:) = 0; u0(end,:) = 0;
bc = sparse(N+2,N+2);

C = 1/(1 + beta*5*pi^2);
analyticalSolution = C.*sin(2*pi* X).*sin(pi*Y);

%% Algebraic solution %%
if 1
    algebraicSolution = zeros(N+2);
    I = speye(N,N);
    subdiag = sparse(2:N,1:N-1,1,N,N);
    D = I/2 - (1/h^2)*beta.*(-2*I + subdiag+subdiag');
    Tdx = kron(D,I)+kron(I,D);
    fhat = reshape(f(2:end-1,2:end-1),N^2,1);
    algebraicSolution(2:end-1,2:end-1) = reshape(Tdx\fhat,N,N);
    algsolerror = norm(algebraicSolution- analyticalSolution)*h;
else
    algsolerror = 1.1800e-06/(4^3);% Define error here the system is big
end

%% FMGV solution %%
range = [6];
timesV = zeros(length(range),1);
timesW = zeros(length(range),1);
errorV = zeros(length(range),1);
errorW = zeros(length(range),1);
iterV = zeros(length(range),1);
iterW = zeros(length(range),1);
itermax = 10;

for ii = 1:length(range)
    disp(['Commencing multigrid iteration with a maximum number of grids: ',num2str(range(ii))])
    MGsol = u0;
    state.gridHistory = [];
    state.numberOfGrids = range(ii);
    state.eta1 = 2;
    state.eta2 = 2;
    state.beta = beta;
    state.nVW = 1;
    currentGrid = 1;
    solutions = cell(1,length(itermax));
    tic
    for jj  = 1:itermax
        disp(['Iteration: ',num2str(jj)])
        [MGsol , state] = FMGV(MGsol , f , currentGrid, state);
        solutions{jj} = MGsol;
    end
    timesV(ii)=+toc/itermax;
    disp(['Average time per multigrid iteration for V-cycle: ', num2str(timesV(ii))])
    disp(['Total time for MG iteration: ', num2str(itermax*timesV(ii))])
    disp('Postprocessing...')
    RMSerror = [norm(u0 - analyticalSolution)*h];
    for jj  = 1:itermax
        disp(jj)
        err = norm(solutions{jj} - analyticalSolution)*h;
        RMSerror = [RMSerror, err];
        if err <= algsolerror
            if iterV(ii) == 0
                iterV(ii) = jj;
            end
        end
    end
    semilogy(0:length(RMSerror)-1,RMSerror);    hold on;
end

T = table(range',timesV, iterV)


%% Comparing with different amount of pre- and postrelaxations %%
