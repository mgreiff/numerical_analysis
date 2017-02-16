function wavesMG(problem)
N = 2^6-1;                         % Number of internal grid points
dt = .01;                          % Time step
h = 1/(N + 1);                     % Spatial step
[X,Y]=meshgrid(linspace(0,1,N+2)); % Mesh of the square [0,1]x[0,1]

switch problem
    case 1 % Initial position as sine function
        u0 = sin(2*pi.*X).*sin(4*pi.*Y);
        v0 = zeros(N+2,N+2);
    case 2 % Initial velocity as sine function
        u0 = zeros(N+2,N+2);
        v0 = 10.*sin(2*pi.*X).*sin(pi.*Y);
    case 3 % Radial wave
        a = 2;
        b = 20;
        u0 = a*exp(-b.*((X-1/2).^2 + (Y-1/2).^2));
        v0 = zeros(N+2,N+2);
    case 4 % Positional pulse in a small area around \approx (0.5,0.5)?0.025
        u0 = zeros(N+2,N+2);
        v0 = zeros(N+2,N+2);
        for ii = 1:N+2
            for jj = 1:N+2
                if sqrt((jj/(N + 2) - 0.5).^2 + (ii/(N + 2)- 0.5).^2) < 0.05
                    v0(ii,jj) = 20;
                end
            end
        end
end

% Enforces boundary conditions
u0(:,1) = 0; u0(:,end) = 0; u0(1,:) = 0; u0(end,:) = 0;

bc = sparse(N+2,N+2);

% Create graphics handle and plot initial state
close all
figure('Name','Waves','Position',[0 0 800 600]);
handle = surf(X,Y,u0);
set(handle,'CDataMapping','direct');

% Set up MG
currentGrid = 1;
state.gridHistory = [];
beta = dt.^2/4;

% Solve problem
uold = u0;
vold = v0;
unew = uold;
vnew = vold;
TdxKernel = [0,1,0;1,-4,1;0,1,0]/h.^2;
TdxPlusKernel = beta.*TdxKernel;
TdxPlusKernel(2,2) = TdxPlusKernel(2,2) + 1;

kineticEnergy = [];

for ii = 1:2000
    % Advance v
    fhat = full(bc);
    fhat(2:end-1,2:end-1) = dt*conv2(uold,TdxKernel,'valid') +...
                               conv2(vold,TdxPlusKernel,'valid');
    MGsol = vnew;
    state.eta1 = 1;
    state.eta2 = 1;
    state.beta = beta;
    state.numberOfGrids = 4;
    state.nVW = 1;

    for jj  = 1:10
        [MGsol , state] = FMGV(MGsol , fhat, currentGrid, state);
    end
    vnew = MGsol;
    % Advance u
    unew(2:end-1,2:end-1) = uold(2:end-1,2:end-1) +...
                            dt*(vnew(2:end-1,2:end-1) +...
                            vold(2:end-1,2:end-1))/2;
    disp(['Iteration: ',num2str(ii),'. Change: ',num2str(norm(unew-uold))])
    vold = vnew;
    uold = unew;
    kineticEnergy = [kineticEnergy, norm(vnew)];
    set(handle,'CData',32*unew+32)
    set(handle,'ZData',unew);
    axis([0 1 0 1 -2 2]);
    drawnow
end
figure(2)
plot(dt.*(1:2000),kineticEnergy/max(kineticEnergy),'b');hold on;
