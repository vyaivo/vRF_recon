% script to run through the reconstruction analysis after the
% fMRI data has been preprocessed and put into MATLAB files.
% these scripts use the full dataset.
% VAV 12/16/2016

sublist = {'AA','AI','AL','AP','AR','AT','AU'};
% voilist = {'V1','V2','V3','V3AB','V4','IPS0'};

bigROIs = {'superVis','superParietal'};
trnFix = 'goodvRFs_TrnAllJitter';

% data was compiled by the following scripts:
% vRFRecon_extractTrialData.m
% vRFRecon_getReconTrialData.m

%%
%%% --- First, do all the stuff for the partial voxel recons (e.g., only
%%% using voxels with well-fit vRFs so we can match to simulations.) --- %%%

% Train & test IEM only using voxels with well-fit vRFs.
% Computes recons & averages across same stim positions.
% trnAllJitter - balanced across L/R runs, as many fix runs as needed
vRFRecon_IEM_goodvRFs_trnAllJitter_acrossSess(sublist,voilist);

% Fit individual subject recons.
tic
for ss = 1:length(sublist)
    thissub = sublist{ss};
    vRFRecon_gridfit_maxRestrict_subRecons(thissub,voilist,trnFix)
end
toc

% Compile all fit data.
compile_recon_data_allsubs(trnFix);

%%
%%% --- Now do the actual simulation part! --- %%%

% Simulate reconstructions under different vRF changes.

% This step takes a long time -- on a very fast computing server, it's ~2500
% seconds in-between file saves for each ROI.
vRFRecon_simulate_vRFchangesToRecon;

%% Make layered IEM data plots.

% Fig 6C: plot RMSE change from null model
% This script will fit all of the recons generated in the previous
% script if you set the doFit flag to 1. Be wary that this step takes even
% longer than the simulation/layered IEM step above -- at least 6 hours on
% a very fast multi-core computing server. If you wish to speed it up, you
% can reduce the number of iterations used to bootstrap CIs on the figures.
simAnalysis_calculateRMSEwFullData;

% Fig 6B: plot amp diff scores for the recon generated with real data (but 
% only good vRFs) vs. diff manipulations of the layered IEM.
simAnalysis_plotModelvsPartialData;

%% --- Control analyses ---
% Now do the simulations /layered IEM with noise scaled by voxel 
% covariance matrix.

% First, calculate voxel covariance based on BOLD data.
simAnalysis_getVoxCovariance_RidgeRegressedResids;

% Now redo the simulation with covaried noise.
vRFRecon_simulate_vRFchangesToRecon_wCovariedNoise;

% Make plot like Fig 6C for the covaried noise
simAnalysis_covModel_calculateRMSEwFullData;






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% redo everything with big ROIs! run this AFTER you have all the individual ROIs already done

% compile big ROI vRFs
compile_vrf_data_allsubs_bigROIs;

vRFRecon_IEM_goodvRFs_trnAllJitter_acrossSess(sublist,bigROIs);

trnFix = 'goodvRFs_TrnAllJitter';
tic
for ss = 1:length(sublist)
    thissub = sublist{ss};
    vRFRecon_gridfit_maxRestrict_subRecons(thissub,bigROIs,trnFix);
end
toc
compile_recon_data_bigROIs(trnFix);

%%% --- Now do the actual simulation part! --- %%%
vRFRecon_simulate_vRFchangesToRecon_bigROIs;
simAnalysis_calculateRMSEwFullData_bigROIs;
simAnalysis_plotModelvsPartialData_bigROIs;