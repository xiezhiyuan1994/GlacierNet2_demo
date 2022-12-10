clear all
close all
clc


% GlacierNet raw output
load('GlacierNet_output')
% raster data 
load('im_test_k2')
% Fused Network
load('Fused_Net')

imt=zeros(size(testing_output));

imt(testing_output>0.5)=1;

imt=bwareaopen(imt,4000);

T1=30;
T2=1000;

imtb=holefill_1(imt,im_c,T1,T2);

% Fused_Network raw output
load('Fused_Net_output')

imt=zeros(size(testing_output));

imt(testing_output>0.5)=1;

imt=bwareaopen(imt,4000);

imtb2=holefill_1(imt,im_c,T1,T2);

imtb_final=imtb2;

imtl2=bwlabel(imtb2);

BW2=double(boundarymask(imtb2));

figure
imshow(imoverlay(uint16(im_c(:,:,[5 4 3])),BW2))  

imt2=terminus_improving(imtb,imtb2,im_c,trainedNet);

figure
imshow(imoverlay(uint16(im_c(:,:,[5 4 3])),boundarymask(imt2)) )


fn2='pakistan-k2-subset-alos-30m.tif';

addpath(genpath('D:\snow_covered_glacier_detect\topotoolbox-master\topotoolbox-master'));

DEM = GRIDobj(fn2);

[DEM,zone] = reproject2utm(DEM,30);

DEM=inpaintnans(DEM,'nearest',1000,8);


DEMf = fillsinks(DEM);

lb7=map_SCG(DEMf,imt2,im_c);
% 
% 
lb7=imresize(lb7,2,'nearest');
% 
% 
bwx=boundarymask(lb7);
% 
figure
imshow(imoverlay(im_c(:,:,[5 4 3]),bwx))