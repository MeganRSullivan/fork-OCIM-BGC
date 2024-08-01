% slurm_run_wrapper.m
% running this script in the slurm matlab call instead of running driver.m
% directly should suppress the command line prompts from appearing in the
% output file
clear all; close all; clc;
diary logs/optPConly_GM15_N23in_npp1_NPPp2c_const.out
%fprintf('----------------------------------- \n\n')
dbstop if error
date_start = datetime('now');
fprintf('Start Date: %s\n',date_start)
fprintf('Optimize the P and C model with oxygen turned off. Goal: to see if the cc and dd parameters give more variable C:P pattern when model is optimized without oxygen. except for O2, SetUp is same as optPCO_GM15_N23in_testNPP. using constant C2P= 117 for satellite NPP conversion. load inital guess on CO and params from hojongs previous run. Temp and DIP_obs are smoothed. \n\n')
%fprintf('continuing the optimization: optimizing the PCO model using the previously optimized cellular growth C:P model from my 2024 GBC paper for both the C2P of uptake, and the p2c conversion for satellite NPP. To see if using consistent C:P in NPP and uptake gives better objective function value. using the SetUp input fields from the Nature paper (Wang et al., 2023) using initial guess on CO and params from from hojongs previous run. no smoothing for Temp obs or DIP_obs \n\n')
run driver_SetUpv2.m
date_end = datetime('now'); fprintf('Complete Date: %s\n',date_end)
fprintf('Total elapsed time was %s\n',date_end - date_start)
diary off