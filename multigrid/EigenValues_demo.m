color = ['g','r','b','k'];
gamma = [1,2/3,0.5,0.25];
beta = 2.5e-5;
N = 2^6-1;
h = 1/(N+1);
figure(1)
for ii = 1:length(gamma)
    f = (1:N)*pi/(N + 1);
    lambda = 1 - gamma(ii) + 2*gamma(ii).*(beta/(h^2+2*beta)).*cos(f);
    subplot(1,2,1); hold on; plot(f,abs(lambda),[color(ii),'o'])
end
legend('\gamma=1','\gamma=2/3','\gamma=1/2','\gamma=1/4')
for ii = 1:length(gamma)
    f = (1:N)*pi/(N + 1);
    lambda = 1 - gamma(ii) + 2*gamma(ii).*(beta/(h^2+2*beta)).*cos(f);
    subplot(1,2,1); hold on; plot(f,abs(lambda),color(ii))
end

title('Spectral radius of the iteration matrix P_{\gamma} with N = 2^6-1 and \beta=1')
xlabel('Normalized frequency')
ylabel('p[P]')

gammasec = [];
factors = 3:15;
Nvec = [];
for ii = factors
    N = 2^ii-1;
    h = 1/(N+1);
    gamma = 1./(1 - 2.*(beta/(h^2+2*beta)).*cos(N*pi/(N + 1)));
    Nvec = [Nvec,N];
    gammasec = [gammasec,gamma];
end
figure(2)
hold on;semilogx(Nvec,gammasec,'+b'); hold on;
hold on;semilogx(Nvec,gammasec,'b--')
T = table(Nvec', gammasec')