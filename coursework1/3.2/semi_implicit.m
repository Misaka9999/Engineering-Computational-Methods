% Solves the equation for the linear spring-mass system with semi-implict Euler method

% Define function name
function y = semi_implicit(x0,v0,N,T,k,m)

% Set default value 
    if nargin==0  
        x0=1,v0=0.1,N=100,T=10,k=5,m=0.5
    end
    
% Set necessary parameters
    taxis  = linspace(0, T, N+1); 
    dt     = T/double(N);

%  Algorithm of semi-implicit Euler method
    
    % Set initial value of v and x
    u_sem = zeros(2,N+1); 
    u0 = [x0;v0];
    u_sem(:,1) = u0

    % Set matrix A
    A = [-k/m*dt 1;-k/m 0 ]; 

    % Iteration of v and x
    for i=1:N
        u_sem(:,i+1) = u_sem(:,i) + dt*A*u_sem(:,i);
    end

    % return value
    y = u_sem  

end