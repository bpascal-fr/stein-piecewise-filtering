%% DENOISING OF PIECEWISE LINEAR SIGNAL
% Sparsity-based denoising for signals, multivariate signals and images
% B. Pascal, N. Pustelnik,  ENS de Lyon, France
% May 2020
%--------------------------------


clear all
close all
clc
addpath(genpath('.'));


% Load an example of real piecewise linear stick-slip signal
load k168v1100.mat
data_crop = data(1:1e3);

% Nonlinear filtering with laplacian
lambda = 5;   
[x_nlf,lambda_nlf,crit_nlf] = nonlinear_filtering(data_crop,'1D','laplacian',lambda,2,1,'direct'); % minimization

% Display signals
figure(3); clf
plot(data_crop,'color',[0.5, 0.5, 0.5])
hold on; grid on
plot(x_nlf,'b')
title('Piecewise linear denoising')
legend('data','denoised')

% Display objective function and duality gap along the iterations
figure(4); clf
subplot(121)
semilogy(crit_nlf.objective);
grid on
title('Objective function')
xlabel('number of iterations')
subplot(122)
semilogy(crit_nlf.gap);
grid on
title('Duality gap')
xlabel('number of iterations')


