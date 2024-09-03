% slurm_run_wrapper.m
% running this script in the slurm matlab call instead of running driver.m
% directly should suppress the command line prompts from appearing in the
% output file
clear all; close all; clc;
diary logs/reoptNature_with_dop_GM15_npp1_gamma20pct.out
%fprintf('----------------------------------- \n\n')
dbstop if error
date_start = datetime('now');
fprintf('Start Date: %s\n',date_start)
fprintf('re-optimize the BGC model from the Nature paper (Wang et al., 2023), but include DOP in the objective function, and scale both DOC and DOP equally. Turn on gamma in SetPar. so gamma = 0.2; fraction of DOP production that goes to the labile pool. Goal: To see how C:P varies as a funciton ofsequestration time if we allow for both labile P and C. the previous version only had labile C, but no labile P, so there was a lot of fast remineralization of carbon, making C:P decrease as a funciton of sequestration time. but this goes against the idea of preferential remineralization of phosporus. I want to see if our result still holds when we have a significant amount of labile DOP. SetUp is same as reoptNature_GM15_npp1_CbPM. using CbPM NPP and GM15 p2c for satellite NPP conversion. Temp and DIP_obs are smoothed. load initial parameters from = ../output/reoptNature_with_dop_GM15_npp1_CTL_He_PCO_DIP1e+00_DIC1e+00_DOC1e+00_ALK1e+00_O21e+00_xhat.mat \n\n')
%fprintf('optimizing the PCO model using the cellular growth C:P model. Use the SetUp input fields from the Nature paper (Wang et al., 2023); use DOP in objective function (to constrain sigP, which affects C:P of POM vs DOM). using initial guess on CO and params from previous reoptNature_with_dop_npp1 run. smoothing both Temp obs and DIP_obs; using CAFE NPP and GM15 p2c for satellite NPP conversion. GOAL: to see if model can find reasonable Cell model C2P and fit to obs when using CAFE NPP compared to CbPM \n\n')
run driver_reoptNature.m
date_end = datetime('now'); fprintf('Complete Date: %s\n',date_end)
fprintf('Total elapsed time was %s\n',date_end - date_start)
diary off