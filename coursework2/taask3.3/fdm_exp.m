function fdm_exp(N,tend,v0,nu,u_0)
% Set default value 
if nargin==0  
    N=64,tend=2.0,v0 = 1.5,nu=0.005, u_0 = @(x) exp(-4*(x-pi).^2);
end


nsteps = 4*N;

L    = 2*pi;
xmesh = linspace(0, L, N+1);
xmesh = xmesh(1:end-1);
dx    = xmesh(2) - xmesh(1);


% create finite difference matrix
e = ones(N,1);
D = spdiags([-e 0*e e], -1:1, N, N);
D(1,N) = -1.0;
D(N,1) = 1.0;

D = 1.0/(2.0*dx)*D;

% define initial value
u0 = u_0(xmesh);


% define the RHS. We use Matlab's fft and ifft function to transfer back
% and forth between physical and spectral space.
% Note how we resolve the nonlinearity by first computing the derivative in
% spectral space through D_x*u, transform u and u_x into physical space using ifft,
% computing u.*u_x and then transforming back into spectral space using
% fft. The linear diffusion term is computed completely in spectral space
% and requires no fft/ifft.

f = @(t,u) (-v0*D*u+nu*D*D*u);
taxis = linspace(0, tend, nsteps+1);
u = exp_euler(u0, tend, nsteps, f);




figure(1);
for n=1:length(taxis)
    
    % plot solution in physical space by applying ifft and taking the real
    % part - note how the solution becomes very steep for small viscosity
    % until a discontinuity in the first derivative forms (a "shock")
    figure(1);
    clf;
    plot(xmesh, u(:,n), 'b'); hold on;
    plot(xmesh, u_0(xmesh), 'r');
    xlim([xmesh(1), xmesh(end)]);
    ylim([-0.2, 1.2]);
    xlabel('x');
    ylabel('u');
    title('\nu=0.005 FDM Solution in physical space');
    txt = strcat('t=',num2str(taxis(n), '%3.2f'));
    text(0.75*L, 0.75, txt);
    drawnow;

    
end

