% slurm_run_wrapper.m
% running this script in the slurm matlab call instead of running driver.m
% directly should suppress the command line prompts from appearing in the
% output file
clear all; close all; clc;
diary logs/reoptNature_Nweighted_dopdoc_GM15_npp1.out
%fprintf('----------------------------------- \n\n')
dbstop if error
date_start = datetime('now');
fprintf('Start Date: %s\n',date_start)
fprintf('re-optimize the BGC model from the Nature paper (Wang et al., 2023), but include DOP in the objective function, and scale both DOC and DOP weights by number of grid cells with obs used in the objective function. Goal: to see if DOM C:P is more reasonable when the model is optimized with less weight on the DOM fields. Except for dopscale and docscale, SetUp is same as reoptNature_GM15_npp1_CbPM. using CbPM NPP and GM15 p2c for satellite NPP conversion. load inital guess on CO and params from hojongs previous run. Temp and DIP_obs are smoothed. \n\n')
%fprintf('continuing the optimization: optimizing the PCO model using the previously optimized cellular growth C:P model from my 2024 GBC paper for both the C2P of uptake, and the p2c conversion for satellite NPP. To see if using consistent C:P in NPP and uptake gives better objective function value. using the SetUp input fields from the Nature paper (Wang et al., 2023) using initial guess on CO and params from from hojongs previous run. no smoothing for Temp obs or DIP_obs \n\n')
run driver_reoptNature.m
date_end = datetime('now'); fprintf('Complete Date: %s\n',date_end)
fprintf('Total elapsed time was %s\n',date_end - date_start)
diary off