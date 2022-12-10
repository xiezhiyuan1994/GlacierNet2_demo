clear all
close all
clc

%% use image stack

load im_test_k2
Z=double(imresize(im_c(:,:,12),0.5,'nearest'));% DEM

r=3;% window size
cellsize=30;% resolution

[geom]=compute_geom_new(Z,r,cellsize);

[SADI]=compute_SADI_new(Z,r,cellsize);


%% use single DEM and topotoolbox

addpath(genpath('D:\snow_covered_glacier_detect\topotoolbox-master\topotoolbox-master'));

DEM = GRIDobj('pakistan-k2-subset-alos-30m.tif');

[geom]=compute_geom2(DEM,r);
