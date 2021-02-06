clear
clc

load('d_inputs.mat')

N = length(Pdrv);

misc.epsilon = 2E4;
misc.maxIterations = 2000;
misc.Emax = Emax;

[E, Pb, time, iters] = f_ADMM(coeffs,Pdrv,Estart,Pbmin,Pbmax,Elowerlim,Eupperlim,P,C,R,V,misc);
fprintf('Time taken using ADMM = %.2f s\n', time)
figure(1)
plot(E)
xlabel('Time (s)')
ylabel('State of Charge (%)')
grid on
legend('ADMM')

try
    [E, Pb, time, iters] = f_CVX(coeffs,Pdrv,Estart,Pbmin,Pbmax,Elowerlim,Eupperlim,P,C,R,V,misc);
    hold on
    plot(E, '--')
    plot([0 N], [Elowerlim Elowerlim], 'r')
    hold off
    legend('ADMM', 'CVX')
catch
    disp('CVX not available')
end



%[dotE, dotmf, Fpred] = f_postprocess(omega, Pdrv, Pb, coeffs, V, R, P, C, engineModel, motorModel);

    

    


