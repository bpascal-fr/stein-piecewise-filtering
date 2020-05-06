%% DENOISING OF PIECEWISE CONSTANT (MULTIVARIATE) SIGNAL
% Sparsity-based denoising for signals, multivariate signals and images
% B. Pascal, N. Pustelnik,  ENS de Lyon, France
% May 2020
%--------------------------------


clear all
close all
clc
addpath(genpath('.'));

% Built a synthetic piecewise constant (multivariate) signal with random
% hops and addtive i.i.d. Gaussian noise
sigma = [1e-1,5e-2,3e-1]; % standard deviation of the noise (optional)
Nval  = [10, 5, 2];       % expected number of differents values taken by the signal

% Generate multivariate data
[data_1D, truth_1D] = create_piecewise_cst_1D_multivariate(3,1000,10,sigma,[10,5,2]);

% Nonlinear filtering with gradient
param.tol   = 1e-3;       % tolerance on duality gap to stop the minimization     (optional)
lambda      = [10,5,30];  % regulariation parameter (global or one per component)
[x_nlf_1D,lambda_nlf_1D,crit_nlf_1D] = nonlinear_filtering(data_1D,'1D','gradient',lambda,2,12,'fourier', param); 

% Display signals
Nl = size(data_1D,1);
figure(1); clf
for nl = 1:Nl
    subplot(1,Nl,nl)
    plot(truth_1D(nl,:),'k')
    hold on; grid on
    plot(data_1D(nl,:),'color',[0.5, 0.5, 0.5])
    plot(x_nlf_1D(nl,:),'b')
    title({'Piecewise constant denoising'})
    xlabel(['Component ',num2str(nl)'])
    legend('truth','noisy','denoised signal')
end

% Display objective function and duality gap along the iterations
figure(2); clf
subplot(121)
semilogy(crit_nlf_1D.objective);
grid on
title('Objective function')
xlabel('number of iterations')
subplot(122)
semilogy(crit_nlf_1D.gap);
grid on
title('Duality gap')
xlabel('number of iterations')

