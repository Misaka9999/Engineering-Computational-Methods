% Solves the equation for the linear spring-mass system with semi-implict Euler method

% Set a few parameter
clear
T      = 10.0;         % final time until which we compute
N      = 10000;          % number of time steps
taxis  = linspace(0, T, N+1); % creates an array 
dt     = T/double(N); % compute the length of each step
k      = 5.0;
b      = 0.0;
m      = 0.5;

% Define the matrix A
A = [0 1 ; -k/m -b];
x0 = 1.0;
v0 = 0.1;
u0 = [x0;v0];

% allocate vectors to store solution; note that for the pendulum the vector
% u has two components.
u_exp      = zeros(2,N+1);
u_exp(:,1) = u0;

% simple forward Euler first
for i=1:N
    u_exp(:,i+1) = u_exp(:,i) + dt*A*u_exp(:,i);
end

% now backward Euler
u_imp = zeros(2,N+1);
u_imp(:,1) = u0;

% define matrix M = Id - dt*A in backward Euler step. IMPORTANT: using 1
% will produce a matrix will all entries equal to one, not the identity.
% Need to use 'eye' instead.
M = eye(2) - dt*A;

for i=1:N  
    u_imp(:,i+1) = M\u_imp(:,i);   
end

% semi_implicit Euler Method
u_sem = semi_implicit(x0,v0,N,T,k,m)


% just for comparison, we also solve it again with Matlab's ode45
f = @(t,u) A*u;
[time_ode45, u_ode45] = ode45(f, [0, T], u0);

% fix size of figure
set(gcf,'Units','centimeter',  'Position',[0 0 18 6]);
%  background color white
set(gcf, 'Color', 'White');

figure(1);
%plot the x
plot(taxis, u_exp(1,:), 'r', 'LineWidth', 2); hold on;
plot(taxis, u_imp(1,:), 'b', 'LineWidth', 2);
plot(taxis, u_sem(1,:), 'g', 'LineWidth', 2);

plot(time_ode45, u_ode45(:,1), 'k', 'LineWidth', 2);
legend('Explicit-x','Implicit-x', 'Semi-x','ode45-x', 'Location','SouthWest');
ylim([-4, 4]);
xlim([taxis(1) taxis(end)]);
xlabel('Time','FontSize',11);
ylabel('x', 'FontSize', 11);

figure(2);
%plot the v
plot(taxis, u_exp(2,:), 'r', 'LineWidth', 2); hold on;
plot(taxis, u_imp(2,:), 'b', 'LineWidth', 2);
plot(taxis, u_sem(2,:), 'g', 'LineWidth', 2);

plot(time_ode45, u_ode45(:,2), 'k', 'LineWidth', 2);
legend('Explicit-v','Implicit-v', 'semi-v','ode45-v', 'Location','SouthWest');
ylim([-5, 5]);
xlim([taxis(1) taxis(end)]);
xlabel('Time','FontSize',11);
ylabel('x', 'FontSize', 11);