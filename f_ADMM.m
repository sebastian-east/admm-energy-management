function [ E, u, time, iterations ] = f_ADMM( coeffs,Pdrv,E0,Pbmin,Pbmax,xmin,xmax,P,C,R,V,misc )

N = length(Pbmin);

rho1 = 2.34E-4;
rho2 = 1E-8;

alpha2 = coeffs(:,1);
alpha1 = coeffs(:,2);
alpha0 = coeffs(:,3);
beta2 = coeffs(:,4);
beta1 = coeffs(:,5);
beta0 = coeffs(:,6);

%% Feasibility check

%% Initialization (these operations can be performed offline and are not timed)

u = zeros(N,1);
u(C) = Pbmin(C);
zeta = zeros(N,1);
lambda1 = zeros(N,1);
lambda2 = zeros(N,1);

I = eye(N);
Psi = tril(ones(N,N));
M = (rho1/rho2 * inv(Psi) * (inv(Psi)') + I);
L = chol(M)';
L = sparse(L);
Diff = sparse(inv(Psi));
Difft = sparse(inv(Psi)');

%% Algorithm
tic

iterations = 0;
flag = 1;

while flag
    
    % The u update for k in P is the solution of a convex optimization
    % problem. This is solved using a newton method that is included in the
    % function f_BacktrackingNewtonVector
    u(P) = f_BacktrackingNewtonVector(alpha0(P), alpha1(P), alpha2(P), beta0(P), beta1(P), beta2(P), V, R, Pdrv(P), rho1, zeta(P), lambda1(P), Pbmin(P), Pbmax(P));
    
    % The x update is trivial
    x = E0 - cumsum(zeta) - lambda2;
    x(x > xmax) = xmax;
    x(x < xmin) = xmin;
    
    %hold the current value of zeta for residual calculations
    zetahold = zeta;
    
    % The zeta update is solved using a method detailed in Appendix D of 
    % "Optimal Power Allocation in Battery/Supercapacitor Electric Vehicles 
    % using Convex Optimization", available at
    % https://ieeexplore.ieee.org/document/9193947
    % and
    % https://arxiv.org/abs/2005.03678
    vec = rho1 * (u + lambda1) - rho2 * cumsum(x - E0 + lambda2, 'reverse');
    vec = vec / rho2;
    vec = Difft * vec; 
    vec = Diff * vec; 
    vec = L \ vec;
    zeta = L' \ vec;
    
    % The residual and Lagrange multiplier updates are also trivial
    r = [u - zeta; x + cumsum(zeta) - E0];
    s = [rho1 * (zetahold - zeta); rho2*cumsum(zetahold - zeta)];
    lambda1 = lambda1 + (u - zeta);
    lambda2 = lambda2 + (x + cumsum(zeta) - E0);
    
    iterations = iterations + 1;
    
    % Termination criteria
    if iterations > misc.maxIterations
        flag = 0;
    end
    
    if max(norm(r), norm(s)) < misc.epsilon
        flag = 0;
    end
    
end

time = toc;

E = E0 - cumsum(u);

return
