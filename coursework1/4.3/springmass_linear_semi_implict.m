% Solves the equation for the linear spring-mass system with semi_implicit Euler method

% Set a few parameter
clear;clc
T      = 10.0;         % final time until which we compute
N      = 100;          % number of time steps
taxis  = linspace(0, T, N+1); % creates an array 
dt     = T/double(N); % compute the length of each step
k      = 5.0;
beta   = 3;
m      = 0.5;
x0     = 1.0;
v0     = 0.1;
u0     = [x0;v0];

% semi_implicit Euler method to nonlinear mass-spring system
% u_sem = semi_implicit(x0,v0,N,T,k,m)
u_sem = semi_implicit_nonlinear(x0,v0,N,T,k,m,beta)

% ode45 to nonlinear mass-spring system
f = @(t,u) [ u(2); -k/m*(u(1) + beta*u(1)^3) ];
[time_ode45, u_ode45] = ode45(f, [0, T], u0);
% fix size of figure
set(gcf,'Units','centimeter',  'Position',[0 0 18 6]);
%  background color white
set(gcf, 'Color', 'White');



figure(1);
%plot x of semi-implicit euler
plot(taxis, u_sem(1,:), 'g', 'LineWidth', 2); hold on;
plot(time_ode45, u_ode45(:,1), 'k', 'LineWidth', 2);
legend('Sem-Implicit-x','ode45-x','Location','SouthWest');
ylim([-2.5, 2.5]);
xlim([taxis(1) taxis(end)]);
xlabel('Time','FontSize',11);
ylabel('x', 'FontSize', 11);

figure(2);
%plot x of semi-implicit euler
plot(taxis, u_sem(2,:), 'g', 'LineWidth', 2); hold on;
plot(time_ode45, u_ode45(:,2), 'k', 'LineWidth', 2);
legend('Sem-Implicit-v','ode45-v','Location','SouthWest');
ylim([-6, 6]);
xlim([taxis(1) taxis(end)]);
xlabel('Time','FontSize',11);
ylabel('x', 'FontSize', 11);
