function [ E,Pb,time,iters ] = f_CVX( coeffs,Pdrv,E0,Pbmin,Pbmax,xmin,xmax,P,C,R,V,misc )


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

%% Algorithm



%length of problem
N = length(Pdrv);

%Generate matrices
Phi = ones(N,1);
Psi = tril(ones(N));

%set up while loop
i = 0;
ticmain = tic;

cvx_begin
    cvx_precision low
    variable Pb(N)
    minimize( cost_function(Pb,Pdrv,alpha2,alpha1,alpha0,beta2,beta1,beta0,V,R) )
    subject to
    Pb <= Pbmax
    Pb >= Pbmin
    Phi*E0 - Psi*Pb <= Phi*xmax
    Phi*E0 - Psi*Pb >= Phi*xmin
    %Phi*E0 - cumsum(Pb) <= Phi*xmax
    %Phi*E0 - cumsum(Pb) >= Phi*xmin
cvx_end

iterations = i;
t = toc(ticmain);
E = E0 - cumsum(Pb);

end

function J = cost_function(Pb,Pdrv,alpha2,alpha1,alpha0,beta2,beta1,beta0,V,R)

ginv = -beta1./2./beta2 + sqrt(-R*Pb.^2./beta2/V^2 + (Pb - beta0)./beta2 + beta1.^2./4./beta2.^2);

Peng = Pdrv - ginv;

Pf = alpha2.*pow_pos(Peng,2) + alpha1.*Peng + alpha0;

J = sum(Pf);
end
