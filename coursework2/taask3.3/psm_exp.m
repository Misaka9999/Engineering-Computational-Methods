function psm_exp(N,tend,v0,nu,u_0)
% Set default value 
if nargin==0  
    N=64,tend=2.0,v0 = 1.5,nu=1, u_0 = @(x) exp(-4*(x-pi).^2);
end

nsteps = 4*N;
assert(mod(N,2)==0, 'N must be even');
L    = 2*pi;
xmesh = linspace(0, L, N+1);
xmesh = xmesh(1:end-1);
D = diag(1j*[0:N/2 -N/2+1:-1]);

u0 = fft(u_0(xmesh).');

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
    plot(xmesh, real(ifft(u(:,n))), 'b'); hold on;
    plot(xmesh, u_0(xmesh), 'r');
    xlim([xmesh(1), xmesh(end)]);
    ylim([-0.2, 1.2]);
    xlabel('x');
    ylabel('u');
    title('\nu=1 PSM Solution in physical space');
    txt = strcat('t=',num2str(taxis(n), '%3.2f'));
    text(0.75*L, 0.75, txt);
    

end

