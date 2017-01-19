function simAnalysis_getVoxCovariance_RidgeRegressedResids(sublist,voi_names,varargin)
% loads data from subjs as well as the best ridge lambda params. calculates
% the residuals between the predicted trial betas for each condition

% last edited 07/30/2016 VAV
% cleaned up for OSF 12/19/2016

%% analysis & path params

if nargin < 1
    sublist = {'AA','AI','AL','AP','AR','AT','AU'};
    voi_names = {'V1','V2','V3','V3AB','V4','IPS0'};
elseif nargin < 2
    voi_names = {'V1','V2','V3','V3AB','V4','IPS0'};
end
%%
root = load_root;
conds = {'Left','Right','Fix'};
nc = length(conds);
simdir = 'sims';
trialdir = 'vRFRecon_trialData';

savefn = fullfile(root,simdir,'ridgeRegression_voxelCovResids.mat');

%% load the old file
for s = 1:length(sublist)
    for v = 1:length(voi_names)
        fname = sprintf('%svRFs/%s_threshVRFRidge_LambdaMinBIC_CVCorr_%s.mat',...
            root,sublist{s},voi_names{v});
        if ~exist(fname,'file')
            disp('No file found. Moving on...');
            continue;
        end
        load(fname);
        fprintf('Loaded vRF data for %s, %s\n',sublist{s},voi_names{v});
        
        %% now get the BOLD data
        
        fn = sprintf('%s%s/%s_vRF_acrossSess_wBalanceInds_%s.mat',root,...
            trialdir,sublist{s},voi_names{v});
        load(fn,'bdat','stimDat','runLabel');        
        % File contains data to estimate vRFs & cross-validate them.
        % bdat:     beta weights for each condition, as a cell array
        % stimDat:  structure containing stimulus info (need for DM)
        % runLabel: run numbers for this condition (L, R, fix)
        
        Xc = cell(1,nc);
        for c = 1:nc
            for rr = 1:length(stimDat{c})
                Xc{c} = cat(1,Xc{c},make_IEM_DM(stimDat{c}(rr),chanX,...
                    chanY,FWHM));
            end
        end

        %% calculate residuals       
        
        % LIST OF ALL VOXEL THRESHOLDING STEPS
        % (1) thresh_r2_vox - voxels with > 0.5 mean R^2 in each session/condition
        % (2) thresh_ridge_vox - voxels that have a minimum BIC
        % (3) thresh_cv_vox - voxels with cross-validation > r_cutoff
        % so to index the original VOI masks, you need to use:
%         bestvox = thresh_r2_vox(thresh_ridge_vox(thresh_cv_vox));
        
        for c = 1:nc
            cbetas = bdat{c}(:,thresh_r2_vox);
            cx = Xc{c};
            cruns = runLabel{c};

            scx = cx ./ repmat(max(cx,[],1),size(cx,1),1);
            w = inv(scx'*scx) * (scx'*cbetas(:,thresh_ridge_vox));

            for i = 1:length(goodlambdas)
                wc(c,:,i) = inv((scx'*scx)+goodlambdas(i)*eye(size(scx,2))) * ...
                    (scx'*cbetas(:,thresh_ridge_vox(i)));
                thisx = (inv( (scx'*scx) + goodlambdas(i)* ...
                    eye(size(scx,2)) ) * scx')';
                thisx = thisx ./ repmat(max(thisx,[],1),size(thisx,1),1);
                b_pred(:,i) = thisx*w(:,i);
                b_rpred(:,i) = thisx*squeeze(wc(c,:,i))';
            end

            % residuals of the LSQ to the trial betas
            realbetas = cbetas(:,thresh_ridge_vox);
            full_bpred{c} = realbetas(:,thresh_cv_vox) - b_pred(:,thresh_cv_vox);
            ridge_bpred{c} = realbetas(:,thresh_cv_vox) - b_rpred(:,thresh_cv_vox);
            
            clear realbetas b_pred b_rpred
        end
        
        beta_pred{s,v} = full_bpred;
        ridge_pred{s,v} = ridge_bpred;
        for c = 1:nc
            ridge_cov{s,v,c} = cov(ridge_bpred{c});
        end
        
        clear bdat thresh_r2_vox thresh_ridge_vox thresh_cv_vox
    end     % end of voi loop
end     % end of sub loop

save(savefn,'beta_pred','ridge_pred','ridge_cov','-v7.3');
fprintf('SAVED %s!\n', savefn);
