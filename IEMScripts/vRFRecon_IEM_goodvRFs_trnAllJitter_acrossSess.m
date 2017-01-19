function vRFRecon_IEM_goodvRFs_trnAllJitter_acrossSess(subList,voi_names,varargin)
% Main reconstruction script. Trains the encoding model on one subset of
% data and tests on the other subset. It then multiplies the reconstructed
% channel weights by the original set of channels, yielding a pixelwise
% reconstruction for each trial. These are then averaged across trials with
% like position.

% VAV 12/7/2016

if nargin < 1
    subList = {'AA','AI','AL','AP','AR','AT','AU'};
    voi_names = {'V1','V2','V3','V3AB','V4','IPS0'};
elseif nargin < 2
    voi_names = {'V1','V2','V3','V3AB','V4','IPS0'};
end
%%
root = '/usr/local/serenceslab/vy/vRFRecon_OSF/';
datadir = 'vRFRecon_trialData';
vdir = 'vRFs';

chanX = linspace(-5,5,9);
chan_spacing = chanX(end)-chanX(end-1);
chanY = (chan_spacing*-5)/2:chan_spacing:(chan_spacing*5)/2;

chan_scale = 1.25;
FWHM = chan_scale*chan_spacing;
plotbool = 0;

warning('error', 'MATLAB:nearlySingularMatrix');
warning('error', 'MATLAB:singularMatrix');

load(sprintf('%s%s/all_vRF_diffscores_noOutliers.mat',root,'vRFfits'),'allout');

%% Loop through each VOI/sub pair, compile data, build & test the spatial IEM
for v = 1:length(voi_names)    
    
    for s = 1:length(subList)  % loop over unique subjects
        skipme = 0;
        
        %% load data
        
        fn = sprintf('%s%s/%s_recon_TrnAllJitter_acrossSess_Bilat-%s.mat',...
            root,datadir,subList{s},voi_names{v});
        load(fn);
        
        % This file has the following variables:
        
        % trnb:         beta weights (nTrainingTrials x nVoxels)
        % trn_c_all:    condition labels (L/R/fix) for all training trials
        % trn_stim_all: stim struct for training trials (to make IEM DM)
        % tst_b_all:    beta weights for test trials
        % tst_c_all:    condition labels for test trials
        % tst_pos_all:  stimulus position for test trials
        
        nconds = length(tst_c_all);
        
        % remove outliers with large difference scores
        switch voi_names{v}
            case 'superVis'                
                load(fullfile(root,'vRFfits','AllSubs_CompiledvRFs_bigROIs.mat'),'keepvox');
                keepvox = keepvox{s,1};
            case 'superParietal'
                load(fullfile(root,'vRFfits','AllSubs_CompiledvRFs_bigROIs.mat'),'keepvox');
                keepvox = keepvox{s,2};
            otherwise
                % load good vRF indices
                vrfn = sprintf('%s%s/%s_threshVRFRidge_LambdaMinBIC_CVCorr_%s.mat',...
                    root,vdir,subList{s},voi_names{v});
                if ~exist(vrfn,'file')
                    fprintf('No vRFs for %s, %s! Skipping...\n',subList{s},voi_names{v});
                    continue;
                end
                load(vrfn,'bestvox');
                keepvox = bestvox; clear bestvox;                
                keepvox(allout{s,v}) = [];
        end
        
        % now only use the good voxels
        trnb = trnb(:,keepvox);
        for c = 1:nconds
            tst_b_all{c} = tst_b_all{c}(:,keepvox);
        end
        
        %% make design matrix
        trnX = [];
        for rr = 1:length(trn_stim_all)
            [tmpX,basis_set,fov,res] = make_IEM_DM(trn_stim_all(rr),chanX,chanY,FWHM);
            trnX = cat(1,trnX,tmpX);
            clear tmpX
        end
        
        %% TRAIN
        trnX = trnX/max(trnX(:));   % need to normalize so weights make sense
        w = trnX\trnb;              % trnb = w * trnX, solve for w
        reconw = w;                 % save out this variable
        
        %% TEST
        for a = 1:nconds
            tst = tst_b_all{a};     % ntrials x nvox

            % then compute the response in each channel on each trial...this gets
            % rid of 'voxel' as a unit and replaces with a representation of the
            % response in each underlying feature channel.
            % this is the "inverted" part of the "inverted encoding model" analysis.
            % this is multivariate - we need all voxels within an ROI 
            
            % TODO: FIX THIS
            try
                x = inv(w*w')*w*tst';
            catch err
                fprintf('%s %s: Too few voxels for matrix inversion.\n',...
                    subList{s},voi_names{v});
                skipme = 1;
                break;
            end
            chanResp{a}  =  x';
            
        end

        %% Compute reconstructions & average across like position.
        
        if skipme
            fprintf('Skipping...\n');
        else
            for c = 1:nconds
                origpos = round(tst_pos_all{c},1);      % original coordinates. need to round out small pixel diffs

                % label all 51 positions separately
                if c == 1
                    posbins = ones(2,size(origpos,1));
                end
                locs = unique(origpos,'rows');
                for l = 1:size(locs,1)
                    lrows = ismember(origpos,locs(l,:),'rows');
                    posbins(c,lrows) = l;
                end
                nbins = unique(posbins(c,:));
                maxbin = max(nbins);
                if c == 1
                    % init some storage arrays
                    nr = size(basis_set,1);             % number of pts in recon
                    ntrials = size(chanResp{c},1);
                    recon_data = nan(nconds,ntrials,nr);
                end
                disp(['Computing recons for condition ' num2str(c)]);
                recon_data(c,:,:) = chanResp{c}*basis_set';

                % now loop through the position bins & average across them.            
                disp('Averaging data...');
                for nb = 1:maxbin
                    recon_avg(c,nb,:) = mean(recon_data(c,posbins(c,:) == nb,:));
                    if plotbool
                        bx = linspace(chanX(1),chanX(end),res(1));
                        by = linspace(chanY(1),chanY(end),res(2));
                        imagesc(bx,by,reshape(recon_avg(c,nb,:),res(2),res(1)));
                        hold on;
                        plot(locs(nb,1),locs(nb,2),'w*');  pause;
                    end
                end

            end     % end loop over conditions        

            %% save output
            stim_pos = tst_pos_all;
            fname = sprintf('%srecons/%s_goodvRFs_TrnAllJitter_acrossSess_%s.mat',...
                root,subList{s},voi_names{v});
            fprintf('Saving %s\n',fname);
            save(fname,'keepvox','chanResp','reconw','chanX','chanY','FWHM','stim_pos',...
                'recon_avg','recon_data','basis_set','fov','res','locs','-v7.3');
        end % end skip me loop
        
    end % end sub loop
end % end VOI loop