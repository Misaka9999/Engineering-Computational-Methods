% Solves the equation for the linear spring-mass system with semi-implict Euler method

% Set a few parameter
clear;clc
T      = 10.0;         % final time until which we compute
% N = [1000, 750, 500, 250, 100, 75, 50];
% N = [2000, 1800, 1600, 1400, 1200, 1000, 800, 600, 400, 200] % list of number of time steps
% N = [20000, 18000, 16000, 14000, 12000, 10000, 8000, 6000, 4000, 2000, 1000,500] % list of number of time steps
N = [200000, 180000, 160000, 140000, 120000, 100000, 80000, 60000, 40000, 20000, 10000,5000] % list of number of time steps
k      = 5.0;
m      = 0.5;
x0 = 1.0;
v0 = 0.1;
u0 = [x0;v0];
omega = sqrt(k/m);

err_semi = zeros(1,length(N));
dts     = zeros(1,length(N));



% define the exact solution as a function handle
u_exact = @(t) x0*cos(omega*t)+v0/omega*sin(omega*t);
v_exact = @(t) -x0*omega*sin(omega*t)+v0*cos(omega*t);

for n=1:length(N)
    taxis = linspace(0, T, N(n)+1); % add +1 to account for t=0
    u_sem = semi_implicit(x0,v0,N(n),T,k,m)
    
    % store the time step dt for plotting
    dts(n) = taxis(2) - taxis(1);
    
    % Now compute the errors
    err_semi(1,n) = max(abs(u_sem(1,:) - u_exact(taxis)));
end

p_semi = polyfit(log(dts), log(err_semi), 1);

figure(1);
loglog(dts, err_semi, 'ro', 'markerfacecolor', 'r'); hold on;
% Because we fitted log(err), we need to apply exp to obtain err.
loglog(dts, exp(polyval(p_semi, log(dts))), 'r-');
txt = strcat('Slope p=', num2str(p_semi(1),'%3.2f'));
text(1e-3, 1e-3, txt); % this needs retuning if values in N are changed
xlim([dts(1) dts(end)])
xlabel('\Delta t');
ylabel('Error');
grid on;
legend('Semi-Implicit Euler', 'Linear fit', 'location', 'NorthWest');


