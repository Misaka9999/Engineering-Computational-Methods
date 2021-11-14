% Solves the equation for the linear spring-mass system semi_implicit Euler method

% Set a few parameter
clear;clc
T      = 10.0;         % final time until which we compute
N      = 1000;          % number of time steps
taxis  = linspace(0, T, N+1); % creates an array 
dt     = T/double(N); % compute the length of each step
k      = 5.0;
m      = 0.5;
omega = sqrt(k/m);
x0 = 1.0;
v0 = 0.1;
u0 = [x0;v0];


% define the exact solution as a function handle
u_exact = @(t) x0*cos(omega*t)+v0/omega*sin(omega*t);
v_exact = @(t) -x0*omega*sin(omega*t)+v0*cos(omega*t);

% calculate the modified discrete energy
u_sem = semi_implicit(x0,v0,N,T,k,m)
E_modified = zeros(1,N+1)

for i=1:N+1
    E_modified(i) = 0.5 * k * u_sem(1,i) * u_sem(1,i) + 0.5 * m * u_sem(2,i) * u_sem(2,i) - 0.5 * k * dt * u_sem(1,i) * u_sem(2,i)
end

%calculate the exact energy
E_exact = zeros(1,N+1);
for i=1:length(taxis)
    E_exact(1:i) = 0.5 * k .* u_exact(taxis(i)) .* u_exact(taxis(i)) + 0.5 * m .* v_exact(taxis(i)) .* v_exact(taxis(i));
end

% fix size of figure
set(gcf,'Units','centimeter',  'Position',[0 0 18 6]);
%  background color white
set(gcf, 'Color', 'White');

figure(1);
%plot x and v in a figure
plot(taxis, E_modified(:), 'b', 'LineWidth', 2);hold on;
plot(taxis, E_exact(:), 'g', 'LineWidth', 2)

legend('Sem-Implicit-Energy','Exact-Energy','Location','SouthWest');
ylim([2.45, 2.55]);
xlim([taxis(1) taxis(end)]);
xlabel('Time','FontSize',11);
ylabel('x', 'FontSize', 11);


