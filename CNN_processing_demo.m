clear all
close all
clc


% Raster data 
load im_test_k2
% Fused network
load('Fused_Net')


testing_output=test_net_for2class_v3(im_c,trainedNet,64,512);

% Fused network raw output
save Fused_Net_output testing_output
