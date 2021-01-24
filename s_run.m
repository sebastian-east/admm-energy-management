clear
clc

load('d_inputs.mat')

misc.epsilon = 2E4;
misc.maxIterations = 2000;
misc.Emax = Emax;

[E, Pb, time, iters] = f_ADMM(coeffs,Pdrv,Estart,Pbmin,Pbmax,Elowerlim,Eupperlim,P,C,R,V,misc);
figure(1)
plot(E)

[E, Pb, time, iters] = f_CVX(coeffs,Pdrv,Estart,Pbmin,Pbmax,Elowerlim,Eupperlim,P,C,R,V,misc);

hold on
plot(E)
plot([0 N], [Eupperlim Eupperlim], 'k')
plot([0 N], [Elowerlim Elowerlim], 'k')
hold off

[dotE, dotmf, Fpred] = f_postprocess(omega, Pdrv, Pb, coeffs, V, R, P, C, engineModel, motorModel);

    

    


