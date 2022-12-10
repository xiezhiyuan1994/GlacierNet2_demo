# GlacierNet2_demo

###### Author: Zhiyuan Xie
###### Year: 2022
###### If you need to use these codes, please cite these papers:
###### GlacierNet: A Deep-Learning Approach for Debris-Covered Glacier Mapping, DOI: 10.1109/ACCESS.2020.2991187
###### GlacierNet2: A Hybrid Multi-Model Learning Architecture for Alpine Glacier Mapping, https://doi.org/10.1016/j.jag.2022.102921

## Please download and add the topotoolbox to the path, such as: (addpath(genpath('topotoolbox path'));
###### The topotoolbox can be download from: https://github.com/wschwanghart/topotoolbox

## The demo data include:
###### 1) im_test_k2: raster data(central Karakoram)
###### 2）GlacierNet_output: GlacierNet raw output
###### 3) Fused_Net_output: Fused network raw output
###### 4) Fused_Net: Fused network
###### 5) pakistan-k2-subset-alos-30m.tif: central Karakoram DEM
###### The above data can be downloaded from: https://drive.google.com/drive/folders/1OdOUUN1Kb0trG5QZNsT4gpGa1nA00Gao?usp=sharing
###### 

## main script files includes
###### 1) CNN_processing_demo.m
###### 2) GlacierNet2_terminus_impoving_SCAZ_estimation_demo.m
###### 3) compute_geom_SADI_demo.m
###### The script 1+2 is the entire GlacierNet2 pipeline. Script 3 is the preprocessing to generate the geomorphometric parameters. Notably, the geomorphometric parameters need to be normalized to the same order of magnitude based on the numerical value range.

## The function files include: 
###### 1) test_net_for2class_v3.m  
###### 2) terminus_improving.m
###### 3) relative_region_merge.m
###### 4) map_SCG.m 
###### 5) landsat8_NDWI.m
###### 6) holefill_1.m
###### 7) compute_SADI_new.m
###### 8) compute_geom2.m
###### 9) compute_geom_new.m
