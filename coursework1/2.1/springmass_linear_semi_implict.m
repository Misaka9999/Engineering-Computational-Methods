% Solves the equation for the linear spring-mass system with semi-implict Euler method

% Set a few parameter
clear;clc
T      = 10.0;         % final time until which we compute
N      = 1000;          % number of time steps
taxis  = linspace(0, T, N+1); % creates an array 
dt     = T/double(N); % compute the length of each step
k      = 5.0;
m      = 0.5;
x0 = 1.0;
v0 = 0.1;
u0 = [x0;v0];

u_sem = semi_implicit(x0,v0,N,T,k,m)

% calculate the energy
E=zeros(1,N+1)
for i=1:N+1
    E(i) = 0.5 * k * u_sem(1,i) * u_sem(1,i) + 0.5 * m * u_sem(2,i) * u_sem(2,i)
end

% fix size of figure
set(gcf,'Units','centimeter',  'Position',[0 0 18 6]);
%  background color white
set(gcf, 'Color', 'White');

figure(1);
%plot x and v in a figure
plot(taxis, E(:), 'b', 'LineWidth', 2);hold on;


legend('Sem-Implicit-Energy','Location','SouthWest');
ylim([2, 3]);
xlim([taxis(1) taxis(end)]);
xlabel('Time','FontSize',11);
ylabel('x', 'FontSize', 11);


