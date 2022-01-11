clear all;
close all;

% basic parameters
tend = 2.0;
v0   = 1.5;
nu   = 0.1;

% Set number of Fourier modes and generate mesh; as always, throw away last
% point which is identical to the first one because of the assumed periodic
% boundary conditions
N = 64;
assert(mod(N,2)==0, 'N must be even');
L    = 2*pi;
xmesh = linspace(0, L, N+1);
xmesh = xmesh(1:end-1);

% create spectral derivative matrix 
D = diag(1j*[0:N/2 -N/2+1:-1]);

% define initial value
% u_0  = @(x) cos(7*x);
u_0 = @(x) cos(3*x);

% compute u by evaluating u0 on the mesh and transforming into spectral
% space
u0 = fft(u_0(xmesh).');
f = @(t,u) (-v0*D*u+nu*D*D*u);

% use ode45 as solver
sol = ode45(f, [0 tend], u0);

fprintf('Minimum time step used: %5.3e \n', min(diff(sol.x)));

% for plotting, interpolate solution to a time mesh we define
taxis = linspace(0, tend, 200);
u = deval(taxis, sol);

% plot the figures of t=2 physical solution and spectral space
figure(1);
subplot(211),plot(xmesh, real(ifft(u(:,length(taxis)))), 'b'); hold on;
plot(xmesh, real(ifft(u(:,1))), 'r');
xlim([xmesh(1), xmesh(end)]);
ylim([-1.1, 1.1]);
xlabel('x');
ylabel('u');
title('Cos(3x) Solution in physical space t = 2');
txt = strcat('t=',num2str(taxis(length(taxis)), '%3.2f'));
text(0.75*L, 0.75, txt);

subplot(212),plot(-N/2:N/2-1, abs(fftshift(u(:,1)))/N, 'rd', 'markerfacecolor', 'r');
plot(-N/2:N/2-1, abs(fftshift(u(:,length(taxis))))/N, 'bo', 'markerfacecolor', 'b');
xlim([-N/2, N/2-1]);
ylim([0.0, 0.55]);
xlabel('k');
ylabel('u_k');
title('Cos(3x) Solution in spectral space t = 2');

