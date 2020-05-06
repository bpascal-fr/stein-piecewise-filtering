%% PIECEWISE CONSTANT IMAGE DENOISING
% Sparsity-based denoising for signals, multivariate signals and images
% B. Pascal, N. Pustelnik,  ENS de Lyon, France
% May 2020
%--------------------------------


clear all
close all
clc

addpath(genpath('.'));


% Load example image
X = double(rgb2gray(imread('cameraman.jpg')));

% Add i.i.d Gaussian noise
Y = X + 10*randn(size(X));

% Piecewise constant denoising
lambda=100;
[x_nlf_2D,lambda_nlf_2D,crit_nlf_2D] = nonlinear_filtering(Y,'2D','gradient',lambda,2,12,'direct'); 

% Display images
figure(5); clf; colormap(pink)
subplot(131)
imagesc(X); axis off image
title('truth')
subplot(132)
imagesc(Y); axis off image
title('noisy')
subplot(133)
imagesc(x_nlf_2D); axis off image
title('denoised')

% Display objective function and duality gap along the iterations
figure(6); clf
subplot(121)
semilogy(crit_nlf_2D.objective);
grid on
title('Objective function')
xlabel('number of iterations')
subplot(122)
semilogy(crit_nlf_2D.gap);
grid on
title('Duality gap')
xlabel('number of iterations')
