%%% vRF_batch_script
%%% calls all vRF analysis functions AND the spatial discrim analysis.

% VAV 3/19/2015
% edited for usability / OSF 12/9/2016

% sub/session & ROI info
sublist = {'AA','AI','AL','AP','AR','AT','AU'};
voilist = {'V1','V2','V3','V3AB','V4','IPS0'};

% data was compiled by the following scripts:
% vRFRecon_extractTrialData.m
% vRFRecon_getvRFTrialData.m
%% Compute vRFs
disp('Calculating ridge parameters...');
vRFRecon_ridgevRFs_wCVThresh(sublist,voilist);

disp('Computing vRFs...');
vRFRecon_computeVRFs(sublist,voilist);

%% Fit vRFs with a surface
disp('Fitting vRFs...');
vRFRecon_gridfit_threshvRFs(sublist,voilist);

%% Compile all data across subjects.
disp('Compiling all vRF data');
% This makes a file with all the vRF fits across all subs/ROIs. It also
% computes the numbers for Table 1.
compile_vrf_data_allsubs;

% Calculate vRF attentional modulations (attend peripheral - attend fix)
vRFanalysis_calcDiffScores_noOutliers;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Make figures/tables & do bootstrapping/randomization statistics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%

% Table 1: vRF threshholding data

root = load_root;
fitdir = 'vRFfits';
load(fullfile(root,fitdir,'AllSubs_CompiledvRFs.mat'));
load(fullfile(root,fitdir,'all_vRF_diffscores_noOutliers.mat'),'allout');

voi_rmse = cell(length(voilist),1);
for v = 1:length(voilist)
    
    noutvox = 0;
    for s = 1:length(sublist)
        % load fit stuff
        fn = sprintf(...
            '%s%s/%s_%s_GridFit_ThreshVRFs_LambdaMinBIC_ALLSESS.mat',...
            root,fitdir,sublist{s},voilist{v});
        if ~exist(fn,'file')
            continue;
        end
        load(fn,'bferr');
        bferr;
        if ~isempty(allout{s,v})
            bferr{3}(allout{s,v},:) = [];
            noutvox = noutvox + numel(unique(allout{s,v}));
        end
        voi_rmse{v} = cat(1,voi_rmse{v},bferr{3});
        clear bferr
    end
    
    vrf_threshhold_table.rmse(v) = mean(voi_rmse{v});
    
    vrf_threshhold_table.thresh4(v) = vrf_threshhold_table.thresh3(v) - ...
        noutvox;
    vrf_threshhold_table.percThresh(v) = vrf_threshhold_table.thresh4(v)...
        ./ vrf_threshhold_table.nvox(v);
end

vrf_threshhold_table.nvox(length(voilist)+1) = sum(vrf_threshhold_table.nvox);
vrf_threshhold_table.thresh1(length(voilist)+1) = sum(vrf_threshhold_table.thresh1);
vrf_threshhold_table.thresh2(length(voilist)+1) = sum(vrf_threshhold_table.thresh2);
vrf_threshhold_table.thresh3(length(voilist)+1) = sum(vrf_threshhold_table.thresh3);
vrf_threshhold_table.thresh4(length(voilist)+1) = sum(vrf_threshhold_table.thresh4);
vrf_threshhold_table.percThresh(length(voilist)+1) = ...
    vrf_threshhold_table.thresh4((length(voilist)+1)) ./ ...
    vrf_threshhold_table.nvox((length(voilist)+1));
vrf_threshhold_table.rmse(length(voilist)+1) = mean(vertcat(voi_rmse{:}));

%%
% Figure 2A: example vRF
vRFanalysis_plotExamplevRF;

% Figure 2B: size/ecc vRF plot
vRFAnalysis_fit_sizeecc_line;

% Figure 2C: vRF vector plot (sample for V4, avg across subs)
% Figure S1: vRF vector plots for attend L/R, all ROIs
vRFAnalysis_makeVRFVectorPlot;

% Figure 2D: mean vRF changes in each ROI
% Figure 2E: polynomial fit lines to all vRF properties
vRFanalysis_attnModulations_byHemisphere_noOutliers_stats;

% Figure S3: vRF coverage plot
vRFanalysis_vRFCoveragePlots;

% Figure S6: attn hemifield x voxel hemisphere
vRFanalysis_attnModulations_byHemisphere_noOutliers_stats;

%%
% Figure 3A: spatial discrim for attend/ignore
% Figure 3B: spatial discrim for p/s/a/b manipulations
vRFRecon_spatialDiscrim;