% script to run through the reconstruction analysis after the
% fMRI data has been preprocessed and put into MATLAB files.
% these scripts use the full dataset.
% VAV 12/01/2016

sublist = {'AA','AI','AL','AP','AR','AT','AU'};
voilist = {'V1','V2','V3','V3AB','V4','IPS0'};

bigROIs = {'superVis','superParietal'};

% data was compiled by the following scripts:
% vRFRecon_extractTrialData.m
% vRFRecon_getReconTrialData.m

%% Train & test IEM. Computes recons & averages across same stim positions.
% trnAllJitter - balanced across L/R runs, as many fix runs as needed
vRFRecon_IEM_trnAllJitter_acrossSess(sublist,voilist);

%% Fit individual subject recons.
tic
for ss = 1:length(sublist)
    thissub = sublist{ss};
    vRFRecon_gridfit_maxRestrict_subRecons(thissub,voilist);
end
toc

%% Compile all fit data.
compile_recon_data_allsubs;

%% Make figures.

% Figure S4 - all averaged recons
% Figure 4B - average recon for AI, V1
% Figure 4C - size/amp plot
reconAnalysis_makePlots;

% Figure 5 - recon parameters for attend/ignore; ANOVA statistics
reconAnalysis_attnVsIgnoreStats;

%% redo everything with super ROIs
vRFRecon_IEM_trnAllJitter_acrossSess(sublist,bigROIs);

tic
for ss = 1:length(sublist)
    thissub = sublist{ss};
    vRFRecon_gridfit_maxRestrict_subRecons(thissub,bigROIs)
end
toc
compile_recon_data_bigROIs;
reconAnalysis_bigROIs_attnVsIgnoreStats;