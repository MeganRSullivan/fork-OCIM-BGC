% slurm_run_wrapper.m
% running this script in the slurm matlab call instead of running driver.m
% directly should suppress the command line prompts from appearing in the
% output file
clear all; close all; clc;
diary logs/transient_test_PC_Atm.out
%fprintf('----------------------------------- \n\n')
dbstop if error
date_start = datetime('now');
fprintf('Start Date: %s\n',date_start)
fprintf('OPTIM = OFF. Run the re-optimized the BGC model from the Nature paper (Wang et al., 2023) that include DOP in the objective function, and scaled both DOC and DOP equally. Set fraction of DOP production that goes to the labile pool, gamma=0. SetUp is same as reoptNature_GM15_npp1_CbPM. using CbPM NPP and GM15 p2c for satellite NPP conversion. Temp and DIP_obs are smoothed. \n\n')
%fprintf('optimizing the PCO model using the cellular growth C:P model. Use the SetUp input fields from the Nature paper (Wang et al., 2023); use DOP in objective function (to constrain sigP, which affects C:P of POM vs DOM). using initial guess on CO and params from previous reoptNature_with_dop_npp1 run. smoothing both Temp obs and DIP_obs; using CAFE NPP and GM15 p2c for satellite NPP conversion. GOAL: to see if model can find reasonable Cell model C2P and fit to obs when using CAFE NPP compared to CbPM \n\n')
run driver_transient_PCcycle.m
date_end = datetime('now'); fprintf('Complete Date: %s\n',date_end)
fprintf('Total elapsed time was %s\n',date_end - date_start)
diary off