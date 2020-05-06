% Sparsity-based denoising for signals, multivariate signals and images
% Stein Unbiased Risk Estimate (SURE) grid search
% Automated tuning of hyperparameters using Stein Unbiased GrAdient
% estimate of the Risk (SUGAR)

% April 2020
% B. Pascal, ENS de Lyon, France
% N. Pustelnik, ENS de Lyon, France
%--------------------------------

clear all
close all
clc

addpath(genpath('.'));


%% 1D UNIVARIATE REAL DATA

load k168v1100.mat
data_1D = data(1:1e3);

param_1D = struct;
param_1D.sigma = std(data_1D(:));

% Grid search
[x_nlf_grid_1D,lambda_grid_1D,sure_vals_1D,lambda_vals_1D,x_grid_1D,lambda_ref_1D, fdmc_1D] = grid_sure(data_1D,'1D','laplacian',-1,1,'direct',param_1D,'log',10);

% BFGS automated tuning
param_1D.fdmc = fdmc_1D; %/!\ it is important for comparison purpose to use the same FDMC parameters for grid search and automated tuning (BFGS)
[x_nlf_1D,lambda_nlf_1D,crit_nlf_1D] = nonlinear_filtering(data_1D,'1D','laplacian',-1,2,1,'direct', param_1D);

figure(1); clf
subplot(121)
plot(data_1D,'color',[0.65,0.65,0.65])
hold on
grid on
plot(x_nlf_grid_1D,'b')
plot(x_nlf_1D,'r')
legend('Signal','Grid search','BFGS (automated)')
title('Estimates')
i_min_sure = find(sure_vals_1D == min(sure_vals_1D));
sure_nlf1_1D = crit_nlf_1D.sure(end);
subplot(122)
semilogx(lambda_vals_1D,sure_vals_1D,'k')
hold on
grid on
semilogx(lambda_vals_1D(i_min_sure),sure_vals_1D(i_min_sure),'k+','markersize',10)
semilogx(lambda_nlf_1D,sure_nlf1_1D,'r*','markersize',10)
xlabel('Regularization parameter \lambda')
legend('SURE','Min. SURE','BFGS (automated)')
title('Grid search v.s. automated tuning')

%% 1D MULTIVARIATE REAL DATA

param_1Dm = struct;
param_1Dm.sigma = [1e-1,5e-2]; % sigma : knwon standard deviation of the noise                          (optional)
datam = create_piecewise_cst_1D_multivariate(2,1e3,10,[1e-1,5e-2],[10,5]);

% Grid search
% One global parameter
% [x_nlf_grid_1Dm,lambda_grid_1Dm,sure_vals_1Dm,lambda_vals_1Dm,x_grid_1Dm,lambda_ref_1Dm, fdmc_1Dm] = grid_sure(datam,'1D','laplacian',-1,12,'direct',param_1Dm,'log',10);
% One parameter per signal component
% /!\ grid search not tractable for more than 2 regularization parameters
[x_nlf_grid_1Dm,lambda_grid_1Dm,sure_vals_1Dm,lambda_vals_1Dm,x_grid_1Dm,lambda_ref_1Dm, fdmc_1Dm] = grid_sure(datam,'1D','laplacian',-2,12,'direct',param_1Dm,'log',10);

% BFGS automated tuning
param_1Dm.fdmc = fdmc_1Dm;
[x_nlf_1Dm,lambda_nlf_1Dm,crit_nlf_1Dm] = nonlinear_filtering(datam,'1D','gradient',-2,2,12,'direct', param_1Dm);

figure(2); clf
subplot(131)
plot(datam(1,:),'color',[0.65,0.65,0.65])
hold on
grid on
plot(x_nlf_grid_1Dm(1,:),'b')
plot(x_nlf_1Dm(1,:),'r')
legend('Signal','Grid search','BFGS (automated)')
title('First component')
subplot(132)
plot(datam(2,:),'color',[0.65,0.65,0.65])
hold on
grid on
plot(x_nlf_grid_1Dm(2,:),'b')
plot(x_nlf_1Dm(2,:),'r')
legend('Signal','Grid search','BFGS (automated)')
title('Second component')
subplot(133); colormap('pink')
imagesc(log(lambda_vals_1Dm(2,:)),log(lambda_vals_1Dm(1,:)),sure_vals_1Dm)
hold on
plot(log(lambda_grid_1Dm(2)),log(lambda_grid_1Dm(1)),'g*')
plot(log(lambda_nlf_1Dm(2)),log(lambda_nlf_1Dm(1)),'gs')
xlabel('\lambda_2')
ylabel('\lambda_1')
title('SURE')
set(gca,'Ydir','normal')
legend('Min. SURE','Automated')


%% 2D REAL IMAGE

X = double(rgb2gray(imread('cameraman.jpg')));
Y = X + 50*randn(size(X));

param_2D = struct;
param_2D.sigma = 50;

% Grid search
[X_nlf_grid_2D,lambda_grid_2D,sure_vals_2D,lambda_vals_2D,X_grid_2D,lambda_ref_2D,fdmc_2D] = grid_sure(Y,'2D','gradient',-1,12,'direct',param_2D,'log',3);

% BFGS automated tuning
param_2D.fdmc = fdmc_2D;
[X_nlf_2D,lambda_nlf_2D,crit_nlf_2D] = nonlinear_filtering(Y,'2D','gradient',-1,2,12,'direct', param_2D);

figure(3);
subplot(131)
imagesc(X); axis off image
title('Noisy image')
subplot(132)
imagesc(X_nlf_grid_2D); axis off image
title('Grid search')
subplot(133)
imagesc(X_nlf_grid_2D); axis off image
title('BFGS (automated)')

figure(6); clf
semilogx(lambda_vals_2D,sure_vals_2D); 
hold on
grid on
semilogy(lambda_nlf_2D,crit_nlf_2D.sure(end),'r+')

