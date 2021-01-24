function [ E,Pb,time,iters ] = f_ADMM( coeffs,Pdrv,E0,Pbmin,Pbmax,xmin,xmax,P,C,R,V,misc )


N = length(Pbmin);

rho1 = 2.34E-4;
rho2 = 8.86E-9;

rho1 = 2.34E-4;
rho2 = 1E-8;

normr = [];
normr2 = [];
norms = [];

alpha2 = coeffs(:,1);
alpha1 = coeffs(:,2);
alpha0 = coeffs(:,3);
beta2 = coeffs(:,4);
beta1 = coeffs(:,5);
beta0 = coeffs(:,6);

%% Feasibility check

%% Algorithm

u = zeros(N,1);
x = zeros(N,1);
zeta = zeros(N,1);
zetahold = zeros(N,1);
lambda1 = zeros(N,1);
lambda2 = zeros(N,1);
r = zeros(2*N,1);
s = zeros(2*N,1);

u(C) = Pbmin(C);

I = eye(N);
Psi = tril(ones(N,N));
M = inv(rho1*I + rho2*(Psi')*Psi);
M2 = (rho1/rho2 * inv(Psi) * (inv(Psi)') + I);
L = chol(M2)';
LT = sparse(L');
L = sparse(L);
Diff = sparse(inv(Psi));
Difft = sparse(inv(Psi)');

PbmaxP = Pbmax(P);
PbminP = Pbmin(P);
alpha2P = alpha2(P);
alpha1P = alpha1(P);
PdrvP = Pdrv(P);

tic

iterations = 0;
flag = 1;

while flag
    
    u(P) = f_BacktrackingNewtonVector(alpha0(P), alpha1(P), alpha2(P), beta0(P), beta1(P), beta2(P), V, R, Pdrv(P), rho1, zeta(P), lambda1(P), Pbmin(P), Pbmax(P));
    
    x = E0 - cumsum(zeta) - lambda2;
    x(x > xmax) = xmax;
    x(x < xmin) = xmin;
    
    vec = rho1 * (u + lambda1) - rho2 * cumsum(x - E0 + lambda2, 'reverse');
    
    zetahold = zeta;
    
    vec = vec / rho2;
    vec = Difft * vec; 
    vec = Diff * vec; 
    vec = L \ vec;
    zeta = L' \ vec;
    
    r = [u - zeta; x + cumsum(zeta) - E0];
    s = [rho1 * (zetahold - zeta); rho2*cumsum(zetahold - zeta)];
    
    lambda1 = lambda1 + (u - zeta);
    lambda2 = lambda2 + (x + cumsum(zeta) - E0);
    
    iterations = iterations + 1;
    
    if iterations > misc.maxIterations
        flag = 0;
    end
    
    if max(norm(r), norm(s)) < misc.epsilon
        flag = 0;
    end
    
end

time = toc;

Pb = u;
E = E0 - cumsum(u);
iters = iterations;

return
