function vRFRecon_ridgevRFs_wCVThresh(subList,voi_names,varargin)
% performs ridge regression on single-voxel RFs. estimates a single lambda
% parameter across all conditions for a single subject. currently collapses
% across day/scan session, leaving one run out in each condition for
% cross-validation.
% lambda is estimated by calculating the BIC for a series of lambda values.
% all voxels which never converge (e.g., BIC never reaches a minimum) are
% discarded. Then we use the mean BIC (across all voxels/conditions) to
% estimate the ridge-regressed vRF.

% created 2/11/2015 VAV
% last edited 11/24/2015 VAV
% cleaned up 6/16/2016 VAV
% cleaned up again 12/5/2016 VAV

if nargin < 1
    subList = {'AA','AI','AL','AP','AR','AT','AU'};
    voi_names = {'V1','V2','V3','V3AB','V4','IPS0'};
elseif nargin < 2
    voi_names = {'V1','V2','V3','V3AB','V4','IPS0'};
end
%%
root = load_root;
trnRuns = {'JitterLeft','JitterRight','JitterFix',...   % all runs are used to estimate vRF
    'DiscreteLeft','DiscreteRight'};                    % but then uses cross-validation
% hemis = {'LH','RH'};
conds = {'Left','Right','Fix'};
nc = length(conds);
trialdir = 'vRFRecon_trialData';

% analysis & path params
chanX = linspace(-5,5,13);              % centers of channels / basis fcns
chan_spacing = chanX(end)-chanX(end-1);
chanY = (chan_spacing*-6)/2:chan_spacing:(chan_spacing*6)/2;
chan_scale = 1.3;
FWHM = chan_spacing*chan_scale;             % width of each channel / basis fcn
lambdas = 0:0.1:750;                        % possible ridge lambda list

% vRF goodness of fit checks
r2_thresh = 0.5;                    % BOLD GLM R^2
r_cutoff = 0.25;                    % cross-validation cutoff

%% 

for v = 1:length(voi_names)
    for s = 1:length(subList)           % loop over unique subjects
        %%        
        fn = sprintf('%s%s/%s_vRF_acrossSess_wBalanceInds_%s.mat',root,...
            trialdir,subList{s},voi_names{v});
        load(fn,'bdat','stimDat','nvox','R2','hemiLabel','runLabel','shufruns');
        
        % File contains data to estimate vRFs & cross-validate them.
        % bdat:     beta weights for each condition, as a cell array
        % stimDat:  structure containing stimulus info (need for DM)
        % R2:       mean GLM R^2 of each voxel on each mapping run
        % runLabel: run numbers for this condition (L, R, fix)
        % shufruns: indices for runs to grab to balance the dataset across
        %           condition (so we have same # for L, R, fix)
        
        fprintf('Calculating vRFs for SUB %s, %s\n', subList{s},voi_names{v});
        
        ntotalvox = nvox; clear nvox;
        
        %% THRESHOLDING STEP 1
        % voxels must pass R^2 beta GLM threshold (r2 > 0.5) across conditions
        for c = 1:nc
            good_r2_vox{c} = find(mean(R2{c},1) > 0.5);
        end
        thresh_r2_vox = intersect(good_r2_vox{1},intersect(good_r2_vox{2},good_r2_vox{3}));
        clear good_r2_vox
        
        %% MAKE DESIGN MATRIX
        % This multiplies a stimulus mask on every trial x bank of 2d
        % cosine filters in order to get a design matrix.
        Xc = cell(1,nc);
        for c = 1:nc
            for rr = 1:length(stimDat{c})
                Xc{c} = cat(1,Xc{c},make_IEM_DM(stimDat{c}(rr),chanX,...
                    chanY,FWHM));
            end
        end
        
        %% FIND THE BEST LAMBDA PARAMETER FOR EACH VOXEL (ESTIMATING ACROSS CONDITIONS)
        
        % First, balance runs across conditions in case this biases
        % the lambda estimation. Also apply the first threshold.
        allb = []; allX = [];
        for c = 1:nc
            thisind = ismember(runLabel{c},shufruns(c,:));
            allb = cat(1,allb,bdat{c}(thisind,thresh_r2_vox));
            allX = cat(1,allX,Xc{c}(thisind,:));
        end
        % normalize DM for each voxel (so vRF units are normalized for each
        % voxel's responses, rather than across all voxels)
        allX = allX ./ repmat(max(allX,[],1),size(allX,1),1);

        % prep parallel computing resources for the next few things
        if isempty(gcp)
            disp('Starting parpool...');
            parpool(8);
        end
        
        % calculate new beta weights for each possible lambda value
        newb = nan(size(allX,2),size(allb,2),length(lambdas));
        rdfs = nan(size(lambdas));
        parfor ri = 1:length(lambdas)
            lamb = lambdas(ri);
            [newb(:,:,ri), rdfs(ri)] = calcRidgeBetas(allb,allX,lamb);
        end

        % calculate Bayesian Information Criterion (BIC) for each lambda
        nb = size(allb,1);
        ridge_bics = nan(size(allb,2),length(lambdas));
        fprintf('calculating BICs for each voxel\n');
        parfor vi = 1:size(allb,2)
            predsignal = allX * squeeze(newb(:,vi,:));
            voxtrn = allb(:,vi);
            [ridge_bics(vi,:)] = calcBIC(predsignal,voxtrn,rdfs,length(lambdas),nb);
        end
        clear predsignal voxtrn

        %% THRESHOLDING STEP 2
        % lambda/BIC curve must have a minimum. if this is not true,
        % the ridge regression solution did not converge.        
        [~, minbici] = min(ridge_bics,[],2);
        thresh_ridge_vox = find(minbici < length(lambdas));
        goodlambdas = lambdas(minbici(thresh_ridge_vox));

        %% EXHAUSTIVELY CROSS VALIDATE EACH VRF AFTER RIDGE REGRESSION
        if isempty(thresh_ridge_vox)
            fprintf('NO VOXELS HAVE MIN BIC! SKIPPING...');
            return;
        end

        wc = nan(nc,size(allX,2),length(goodlambdas));
        for c = 1:nc
            cbetas = bdat{c}(:,thresh_r2_vox);
            cx = Xc{c};
            cruns = runLabel{c};

            scx = cx ./ repmat(max(cx,[],1),size(cx,1),1);
            for i = 1:length(goodlambdas)
                wc(c,:,i) = inv((scx'*scx)+goodlambdas(i)*eye(size(scx,2))) * ...
                    (scx'*cbetas(:,thresh_ridge_vox(i)));
            end

            %w = inv(scx'*scx) * (scx'*cbetas(:,thresh_ridge_vox));

            nr = length(unique(cruns));
            ridge_r = nan(nc,nr,length(goodlambdas));
            % CALCULATE CROSS-VALIDATION CORRELATION ON LEFT OUT RUN
            fprintf('Cross-validating condition: %s\n', conds{c});
            for lr = 1:nr
                % data from most runs
                b = cbetas(cruns ~= lr,:);
                X = cx(cruns ~= lr,:);
                X = X ./ repmat(max(X,[],1),size(X,1),1);   % scale channels from 0 - 1
                for i = 1:length(goodlambdas)
                    thislambda = goodlambdas(i);
                    w(:,i) = inv((X'*X)+thislambda*eye(size(X,2)))*...
                        (X'*b(:,thresh_ridge_vox(i)));
                end
                % data from left out run 
                lb = cbetas(cruns == lr,:);
                lx = cx(cruns == lr,:);
                lx = lx ./ repmat( max(lx,[],1),size(lx,1),1 );        
                for i = 1:length(goodlambdas)
                    thisx = (inv( (lx'*lx) + goodlambdas(i)* ...
                        eye(size(lx,2)) ) * lx')';
                    thisx = thisx ./ repmat(max(thisx,[],1),size(thisx,1),1);
                    b_ridge(:,i) = thisx*w(:,i);
                    goodlb(:,i) = lb(:,thresh_ridge_vox(i));
                    clear thisx
                end

                % calculate correlation with the full data
                ridge_r(c,lr,:) = diag(corr(goodlb,b_ridge));                
                clear goodlb b_ridge b X w lb lx
            end
            % FIND MEAN CROSS VALIDATION R FOR EACH VOXEL
            meanr(c,:) = nanmean(ridge_r(c,1:nr,:),2);
        end

        %% THRESHOLDING STEP 3
        % RIDGE REGRESSION SOLUTIONS FOR EACH VOXEL MUST PASS A 
        % CROSS-VALIDATION THRESHOLD OF AT LEAST 0.25        
        
        for c = 1:nc
            cutoffvox{c} = find(meanr(c,:) >= r_cutoff);
        end
        thresh_cv_vox = intersect( intersect(cutoffvox{1},cutoffvox{2}), cutoffvox{3} );
        for cc = 1:nc
            w_conds(cc,:,:) = squeeze(wc(cc,:,thresh_cv_vox));
        end

        % LIST OF ALL VOXEL THRESHOLDING STEPS
        % (1) thresh_r2_vox - voxels with > 0.5 mean R^2 in each session/condition
        % (2) thresh_ridge_vox - voxels that have a minimum BIC
        % (3) thresh_cv_vox - voxels with cross-validation > r_cutoff
        % so to index the original VOI masks, you need to use:
        bestvox = thresh_r2_vox(thresh_ridge_vox(thresh_cv_vox));
        
        % save out hemisphere labels for the best voxels        
        bestvoxhemi = hemiLabel(bestvox);

        %% save out data
        fname = sprintf('%svRFs/%s_threshVRFRidge_LambdaMinBIC_CVCorr_%s.mat',...
            root,subList{s},voi_names{v});
        fprintf('saving to %s...\n\n',fname);
        save(fname,'ntotalvox','bestvox','thresh_r2_vox','thresh_ridge_vox',...
            'thresh_cv_vox','w_conds','goodlambdas','r_cutoff','chanX','chanY',...
            'FWHM','lambdas','ridge_r','R2','cutoffvox','hemiLabel',...
            'bestvoxhemi','-v7.3');
        
        clear meanr bestvox thresh_r2_vox thresh_ridge_vox thresh_cv_vox 
        clear w_conds goodlambdas cutoffvox
        
    end         % end of sub loop
end             % end of voi loop

end     % function end


function [bic] = calcBIC(predsignal,voxtrn,df,lsize,tsize)

ssres = sum((repmat(voxtrn,1,lsize) - predsignal).^2, 1);
bic = tsize*log(ssres) + df*log(tsize);

end

function [ww, df] = calcRidgeBetas(y,X,lambda)
% here uses inv(X'X+lambda*eye)*X'y to find ridge betas at all voxels

idmat = eye(size(X,2));
ww = inv( (X'*X)+lambda*idmat )*(X'*y);
df = trace(X*inv((X'*X)+lambda*idmat)*X');

end