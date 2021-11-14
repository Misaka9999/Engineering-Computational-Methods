% Solves the equation for the linear spring-mass system with semi-implict Euler method

% Set a few parameter
clear;clc
T      = 10.0;         % final time until which we compute
N      = 100;          % number of time steps
taxis  = linspace(0, T, N+1); % creates an array 
dt     = T/double(N); % compute the length of each step
k      = 5.0;
m      = 0.5;
x0 = 1.0;
v0 = 0.1;
u0 = [x0;v0];

u_sem = semi_implicit(x0,v0,N,T,k,m)

% fix size of figure
set(gcf,'Units','centimeter',  'Position',[0 0 18 6]);
%  background color white
set(gcf, 'Color', 'White');

figure(1);
%plot x and v in a figure
plot(taxis, u_sem(1,:), 'g', 'LineWidth', 2);hold on;
plot(taxis, u_sem(2,:), 'b', 'LineWidth', 2);

legend('Sem-Implicit-x','Sem-Implicit-v','Location','SouthWest');
ylim([-4, 4]);
xlim([taxis(1) taxis(end)]);
xlabel('Time','FontSize',11);
ylabel('x', 'FontSize', 11);

figure(2);
%plot x of semi-implicit euler
plot(taxis, u_sem(1,:), 'g', 'LineWidth', 2); hold on;

legend('Sem-Implicit-x','Location','SouthWest');
ylim([-2.5, 2.5]);
xlim([taxis(1) taxis(end)]);
xlabel('Time','FontSize',11);
ylabel('x', 'FontSize', 11);

figure(3);
%plot x of semi-implicit euler
plot(taxis, u_sem(2,:), 'g', 'LineWidth', 2); hold on;

legend('Sem-Implicit-v','Location','SouthWest');
ylim([-4, 4]);
xlim([taxis(1) taxis(end)]);
xlabel('Time','FontSize',11);
ylabel('x', 'FontSize', 11);
