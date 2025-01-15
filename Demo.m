%% =================================================================
% This script runs the FCETN decomposition-based TC method
%
% More detail can be found in [1]
% [1] Yao Li, Yujie Zhang, Hongwei Li
%     Hyperspectral Image Completion Using
%   Fully-Connected Extended Tensor Network
%      Decomposition and Total Variation
%
% Please make sure your data is in range [0, 1].
%


%% =================================================================
clc;
clear;
close all;
addpath(genpath('lib'));
addpath(genpath('data'));
addpath('Tools/functions');
addpath('Tools/prox_operators');

%%
methodname    = {'Observed', 'FCETN-TV'};
Mnum          = length(methodname);
Re_tensor     = cell(Mnum,1);
psnr          = zeros(Mnum,1);
ssim          = zeros(Mnum,1);
time          = zeros(Mnum,1);

%% Load initial data
load('Africa.mat')
load('mask.mat')
T=Africa;
%% Observed tensor
Ndim      = ndims(T);
Nway      = size(T);
Omega     =find(mask==1);
F         = zeros(Nway);
F(Omega)  = T(Omega);
%%
i  = 1;
Re_tensor{i} = F;
[psnr(i), ssim(i)] = quality_ybz(T*255, Re_tensor{i}*255);
enList = 1;

%% Perform  algorithms

i = i+1;
% initialization of the parameters
% Please refer to our paper to set the parameters
opts=[];
opts.max_R = [0, 8, 8;
              0, 0, 8;
              0, 0, 0];  
opts.RR    = [0, 1, 1;
              0, 0, 1;
              0, 0, 0];
opts.tol   = 1e-4;
opts.maxit =1000;
opts.rho   = 0.001;
opts.Frame    = 1; % (0,1,3)
opts.Level    = 1;  % [1,2,3,4,5,6]
opts.Lamda=1;
opts.mu=1.1;
opts.beta=1.1;
%%%%%
fprintf('\n');
disp(['performing ',methodname{i}, ' ... ']);
t0= tic;
[Re_tensor{i},G,Out]        = inc_FCETN_TV(F,Omega,opts);
time(i)                     = toc(t0);
[psnr(i), ssim(i)]          = quality_ybz(T*255, Re_tensor{i}*255);
enList = [enList,i];

%% Show result
fprintf('\n');
fprintf('================== Result =====================\n');
fprintf(' %8.8s    %5.4s    %5.4s    \n','method','PSNR', 'SSIM' );
for i = 1:length(enList)
    fprintf(' %8.8s    %5.3f    %5.3f    \n',...
    methodname{enList(i)},psnr(enList(i)), ssim(enList(i)));
end
sam=SAM3D(T*255, Re_tensor{2}*255);

X4=mask(:,:,5);
imshow(X4/max(X4(:)),'border','tight','initialmagnification','fit')
set(gcf,'Position',[0,0,256,256]);
axis normal;