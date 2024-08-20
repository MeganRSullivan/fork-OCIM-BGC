% slurm_run_wrapper.m
% running this script in the slurm matlab call instead of running driver.m
% directly should suppress the command line prompts from appearing in the
% output file
clear all; close all; clc;
diary logs/optPCO_Cell_with_dop_N23in_npp2_NPPp2cGM15.out
%fprintf('----------------------------------- \n\n')
dbstop if error
date_start = datetime('now');
fprintf('Start Date: %s\n',date_start)
%fprintf('re-optimize the BGC model from the Nature paper (Wang et al., 2023), but include DOP in the objective function, and scale both DOC and DOP weights by number of grid cells with obs used in the objective function. Goal: to see if DOM C:P is more reasonable when the model is optimized with less weight on the DOM fields. Except for dopscale and docscale, SetUp is same as reoptNature_GM15_npp1_CbPM. using CbPM NPP and GM15 p2c for satellite NPP conversion. load inital guess on CO and params from hojongs previous run. Temp and DIP_obs are smoothed. \n\n')
fprintf('optimizing the PCO model using the cellular growth C:P model. Use the SetUp input fields from the Nature paper (Wang et al., 2023); use DOP in objective function (to constrain sigP, which affects C:P of POM vs DOM). using initial guess on CO and params from previous reoptNature_with_dop_npp1 run. smoothing both Temp obs and DIP_obs; using CAFE NPP and GM15 p2c for satellite NPP conversion. GOAL: to see if model can find reasonable Cell model C2P and fit to obs when using CAFE NPP compared to CbPM \n\n')
run driver_SetUpv2.m
date_end = datetime('now'); fprintf('Complete Date: %s\n',date_end)
fprintf('Total elapsed time was %s\n',date_end - date_start)
diary off