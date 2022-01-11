function psm_exp(N,tend,v0,nu,u_0)

% Set default value 
if nargin==0  
    N=64,tend=2.0,v0 = 1.5,nu=1, u_0 = @(x) exp(-4*(x-pi).^2);
end
% Set number of Fourier modes and generate mesh; as always, throw away last
% point which is identical to the first one because of the assumed periodic
% boundary conditions
nsteps = 4*N;
assert(mod(N,2)==0, 'N must be even');
L    = 2*pi;
xmesh = linspace(0, L, N+1);
xmesh = xmesh(1:end-1);
D = diag(1j*[0:N/2 -N/2+1:-1]);
% compute u by evaluating u0 on the mesh and transforming into spectral
% space
u0 = fft(u_0(xmesh).');

f = @(t,u) (-v0*D*u+nu*D*D*u);
taxis = linspace(0, tend, nsteps+1);

% use explict euler method as solver
u = exp_euler(u0, tend, nsteps, f);

end