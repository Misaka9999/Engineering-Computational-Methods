function fdm_exp(N,tend,v0,nu,u_0)
% Set default value 
if nargin==0  
    N=64,tend=2.0,v0 = 1.5,nu=0.005, u_0=@(x) cos(x); %u_0 = @(x) exp(-4*(x-pi).^2);
end

% set basic parameters
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

end
