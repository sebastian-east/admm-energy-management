function [ Pb ] = f_BacktrackingNewtonVector(alpha0,alpha1,alpha2,beta0,beta1,beta2,V,R,Pdrv,rho1,zeta,lambda1,Pbmin,Pbmax)
% Dimensions of input arrays (all double unless otherwise specified) -
%
% alpha0 = |P| x 1
% alpha1 = |P| x 1
% alpha1 = |P| x 1
% beta0 = |P| x 1
% beta1 = |P| x 1
% beta2 = |P| x 1
% V = 1 x 1
% R = 1 x 1
% Pdrv = |P| x 1
% rho2 = 1 x 1
% rho3 = 1 x 1
% zeta = |P| x 1
% eta = |P| x 1
% lambda2 = |P| x 1
% lambda3 = |P| x 1
% Pbmax = |P| x 1

%Dimensions of output arrays - 
% Pb = N x 1

%% Initialise decision variable
x = Pbmax;

%Calculate coefficients
d = -beta1./2./beta2;
a = -R./beta2./V.^2;
b = 1./beta2;
c = -beta0./beta2 + beta1.^2./4./beta2.^2;

%Calculate derivatives

% ginv = d + sqrt(a.*x.^2 + b.*x + c); %Explicit derivatives
% ginvdash = 0.5.*(2.*a.*x+b).*(a.*x.^2 + b.*x + c).^-0.5;
% ginvddash = -0.25.*(2.*a.*x+b).^2.*(a.*x.^2 + b.*x + c).^-1.5 + 0.5.*2.*a.*(a.*x.^2 + b.*x + c).^-0.5;

sqt_quad = realsqrt(a.*x.^2 + b.*x + c); %Derivative calculations optimized for speed
ginv = d + sqt_quad;
ginvdash = 0.5.*(2.*a.*x + b)./sqt_quad;
ginvddash = -ginvdash.^2./sqt_quad + a./sqt_quad;

h = Pdrv - ginv;

%f = alpha2.*h.^2 + alpha1.*h + alpha0;
fdash = 2.*alpha2.*h + alpha1;
fddash = 2.*alpha2;

%y = f + rho1/2*(x - zeta + lambda_2).^2;
ydash = -ginvdash.*fdash + rho1.*(x - zeta + lambda1);
yddash = (ginvdash).^2.*fddash - ginvddash.*fdash + rho1;

%Backtracking parameter
beta = 0.5;

%Start Newton method with Newton step
deltax = -ydash./yddash;

iter = 1;
tstart = ones(size(deltax));

flag = 1;

while flag
    
    %Initialise iteration
    t = tstart;
    iter2 = 1;
    
    %Calculate Newton step
    deltax = -ydash./yddash;
    
%     %Check ginv is real
%     ginv2 = a.*(x + t.*deltax).^2 + b.*(x + t.*deltax) + c;
%     backtrack = ginv2 < 0;
%     
%     %Backtrack until ginv is real
%     while max(backtrack) > 0
%         
%         t(backtrack) = beta*t(backtrack);
%         ginv2 = a.*(x + t.*deltax).^2 + b.*(x + t.*deltax) + c;
%         backtrack = ginv2 < 0;
%         
%         iter2 = iter2 + 1
%         
%         if iter2 > 2000
%             disp('broken')
%             return
%         end
%     end
    
    %Take Newton Step
    xhold = x;
    x = x + t.*deltax;
    x = min(Pbmax, max(Pbmin, x));
    ginv2 = a.*(x).^2 + b.*(x) + c;
    
    %Update derivatives
    
    %Explicit derivatives
    %     ginv = d + sqrt(a.*x.^2 + b.*x + c);
    %     ginvdash = 0.5.*(2.*a.*x + b).*(a.*x.^2 + b.*x + c).^-0.5;
    %     ginvddash = -0.25.*(2.*a.*x+b).^2.*(a.*x.^2 + b.*x + c).^-1.5 + 0.5.*2.*a.*(a.*x.^2 + b.*x + c).^-0.5;
    
    %Optimized derivatives for computational speed
    sqt_quad = realsqrt(ginv2);
    ginv = d + sqt_quad;
    ginvdash = 0.5.*(2.*a.*x + b)./sqt_quad;
    ginvddash = -ginvdash.^2./sqt_quad + a./sqt_quad;
    
    h = Pdrv - ginv;
    
    %f = alpha2.*h.^2 + alpha1.*h + alpha0;
    fdash = 2.*alpha2.*h + alpha1;
    fddash = 2.*alpha2;
    
    %y = f + rho1/2*(x - zeta + lambda_2).^2;
    ydash = -ginvdash.*fdash + rho1.*(x - zeta + lambda1);
    yddash = (ginvdash).^2.*fddash - ginvddash.*fdash + rho1;
    
    iter = iter + 1;
    iter2 = iter2 + 1;
    
    if norm(x - xhold) < 1E-5
        flag = 0;
    end
    
end

Pb = x;

end

