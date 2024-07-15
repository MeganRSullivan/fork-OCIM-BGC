% slurm_run_wrapper.m
% running this script in the slurm matlab call instead of running driver.m
% directly should suppress the command line prompts from appearing in the
% output file
clear all; close all; clc;
diary logs/optPCO_Const_prescribe_C2P_N23in_NPPp2c_Const_optGBC2024.out
%fprintf('----------------------------------- \n\n')
dbstop if error
date_start = datetime('now');
fprintf('Start Date: %s\n',date_start)
%fprintf('Reoptimizing Temperature-dependent C2P model using the SetUp input processing from the Nature paper (Wang et al., 2023) (smooth WOA PO4 and Temp obs; remove arctic docraw,alkraw,dicraw ; remove med docraw, alkraw, dicraw, po4raw, o2raw. ; rmoutliers from docraw, po4raw, dicraw, alkraw, o2raw); using WOA18 data files created in make_datafile_24layer/bin_WOA_nut.m (units converted using constant seawater density); using initial guess on CO from from hojongs previous run. initial parameter values from Hojongs previous run. NPP p2c conversion = 1/117 \n\n')
fprintf('optimizing the PCO model using the previously optimized constant C:P model from my 2024 GBC paper for both the C2P of uptake, and the p2c conversion for satellite NPP. To see if using consistent C:P in NPP and uptake gives better objective function value. using the SetUp input fields from the Nature paper (Wang et al., 2023) using initial guess on CO and params from from hojongs previous run. no smoothing for Temp obs or DIP_obs \n\n')
run driver_SetUpv2.m
date_end = datetime('now'); fprintf('Complete Date: %s\n',date_end)
fprintf('Total elapsed time was %s\n',date_end - date_start)
diary off