function vRFRecon_simulate_vRFchangesToRecon_wCovariedNoise(doROI)
% giant script to run through all the possible / desired combos of vrf
% changes. also calculates an error metric.
% compute variance in each recon pixel to normalize the SSE. also just
% saves out raw SSE & RMSE.

% VAV 5/12/2016
% started OSF cleanup 12/14/2016

if nargin < 1
    doROI = 1:6;
end
%%
root = load_root;
fdir = 'vRFfits';
rdir = 'recons';
simdir = 'sims';
behavdir = 'vRFRecon_behavData';
trialdir = 'vRFRecon_trialData';
pdir = 'vRFs';
clist = {'Left','Right','Fix'};
nconds = length(clist);

vrf_w_std = 0.5;
sim_response_std = 0.5;

parcombs = [1 1 1 1 1;
            1 1 0 1 1;
            1 1 1 0 1;
            0 0 1 1 1;
            1 1 0 1 0;
            0 0 1 1 0;
            1 1 1 0 0;
            1 1 0 0 0;
            0 0 0 1 0;
            0 0 1 0 0;
            0 0 0 0 0];
        
nm = size(parcombs,1);
niters = 100;


%% load in & setup some vars/data

% load vRF data
vrfn = fullfile(root,fdir,'AllSubs_CompiledvRFs.mat');
load(vrfn);

% load voxel covariance matrix
load(fullfile(root,simdir,'ridgeRegression_voxelCovResids.mat'),'ridge_cov');

% get vRF outliers (based on vRF diff scores > 3*sd from mean)
voutfn = sprintf('%s%s/all_vRF_diffscores_noOutliers.mat',root,fdir);
load(voutfn,'allout');

% define some params based on the first sub/VOI
exampleFile = sprintf('%s%s/%s_TrnAllJitter_acrossSess_V1.mat',...
    root,rdir,sublist{1});
load(exampleFile,'chanX','chanY','res','fov','FWHM','locs');
npix = res(1)*res(2);

% define the resolution of the simulated vRF activity using first sub
% this should be much smaller than the IEM matrix
stepsize = 0.1;
xlist = min(chanX):stepsize:max(chanX);
ylist = min(chanY):stepsize:max(chanY);    
[rx, ry] = meshgrid(xlist,ylist);
xx = rx(:);
yy = ry(:);
simsize = (1/rad2fwhm(1))*1.3*(xlist(2)-xlist(1));

% also make an encoding model matrix
% basis_size = (1/rad2fwhm(1))*FWHM;
[cx, cy] = meshgrid(chanX,chanY);
basis_size = (1/rad2fwhm(1)) * FWHM;
smallres = [round((length(chanX)/length(chanY)) * 51) 51];
if mod(smallres(1),2) == 0
    smallres(1) = smallres(1)+1;
end
gridxPts = linspace(-fov(1)/2,fov(1)/2,smallres(1));
gridyPts = linspace(-fov(2)/2,fov(2)/2,smallres(2));
[gridx, gridy] = meshgrid(gridxPts,gridyPts);
basis_set = build_basis_pts(cx(:),cy(:),basis_size,gridx(:),gridy(:),0)';

%%

for voin = doROI
    skipsub = [];
    randseed = (sum(clock)*10000);
    rng(randseed);
    
    savefn = sprintf('%s%s/sim_covNoise_vRFRecon_allSubsNoOutliers_%s.mat',root,simdir,...
        voilist{voin});
    tic
    for sn = 1:ns
        sub = sublist{sn};
        subdat = alldat{sn,voin};        
        subdat(allout{sn,voin},:,:) = [];       % remove outliers
        % need to convert vRF FWHM size to size constant for make2dcos fcn
        subdat(:,:,3) = subdat(:,:,3) ./ rad2fwhm(1);
        
        % also remove outliers in covariance matrix
        subcov = {ridge_cov{sn,voin,:}};
        for c = 1:nconds
            subcov{c}(allout{sn,voin},:) = [];
            subcov{c}(:,allout{sn,voin}) = [];
        end

        if isempty(subdat) || size(subdat,1) == 1 || ismember(sn,skipsub)
            % if there's not enough data or we should skip the subject, skip
            if ~ismember(sn,skipsub)
                skipsub = [skipsub sn];
            end
            continue;
        end
        nvox = size(subdat,1);

        tfn = sprintf('%s%s/%s_recon_TrnAllJitter_acrossSess_Bilat-%s.mat',...
            root,trialdir,sub,voilist{voin});
        load(tfn,'trn_stim_all','tst_pos_all');
        if ~exist('npos','var')
            nconds = nc; clear nc;
            npos = size(trn_stim_all(1).locsGridDeg,1);
            subrecon = nan(ns,niters,nm,nconds-1,npos,size(basis_set,2));
            all_sim_recon = nan(niters,nm,nconds-1,npos,size(basis_set,2));
        end

        % copy over some basic stim params before we write in real values
        trnstim = trn_stim_all(1);
        tststim = trnstim;
        uniqPos = locs;

        % generate training dataset params using positions of jittered stimuli
        trn_condList = nan(size(trn_stim_all,1)*npos,1);
        idlast = 1;
        for rr = 1:length(trn_stim_all)
            % now copy all the run-specific stim locations
            idx = idlast:idlast+npos-1;
            trn_condList(idx) = trn_stim_all(rr).conditions;
            trnstim.locsRealDeg(idx,:) = trn_stim_all(rr).locsRealDeg;
            trnstim.xLocDeg(idx) = trn_stim_all(rr).locsRealDeg(:,1);         
            trnstim.yLocDeg(idx) = trn_stim_all(rr).locsRealDeg(:,2);
            idlast = idlast+npos;
        end

        % for the test data, generate a separate list from the discrete stimuli
        nTestTrials = size(tst_pos_all{1},1);
        nTestRuns = nTestTrials ./ npos;
        test_list = nan(nTestTrials*(nconds-1),2);              % stim ID, condition
        idlast = 1;
        for cc = 1:nconds-1
            for rr = 1:nTestRuns
                idx = idlast:idlast+npos-1;
                test_list(idx,1) = 1:npos;
                test_list(idx,2) = repmat(cc,npos,1);
                tststim.xLocDeg(idx) = uniqPos(:,1);         
                tststim.yLocDeg(idx) = uniqPos(:,2);
                idlast = idlast + npos;
            end
        end
        
        fprintf('Making vRF simulation DMs...\n');
        % This is setup to simulate the vRF responses on each trial.
        trn_all = sim_RF_DM(trnstim,xlist,ylist,simsize,fov);
        tst_all = sim_RF_DM(tststim,xlist,ylist,simsize,fov);
        % this took 155 seconds on nc4 with 8 cores

        [trnX,~,~,~] = make_IEM_DM(trnstim,chanX,chanY,FWHM,fov,res);
        trnX = trnX ./ max(trnX(:));

        %% now that data is loaded, iterate through with noise sims        
        iterrecon = nan(niters,nm,nconds-1,npos,size(basis_set,2));
        parfor ii = 1:niters
%         for ii = 1:niters
            disp(ii);
            tmp = doSimRecon(parcombs,subdat,nconds,nvox,vrf_w_std,...
                sim_response_std,trn_condList,test_list,trn_all,tst_all,...
                trnX,basis_set,xx,yy,subcov);
            if isempty(tmp)
                skipsub = [skipsub,sn];
                continue;
            end
            iterrecon(ii,:,:,:,:) = tmp;
        end
        subrecon(sn,:,:,:,:,:) = iterrecon;
        
        clear subdat subcov
    end
    toc
    all_sim_recon = squeeze(nanmean(subrecon,1));
    save(savefn,'all_sim_recon','skipsub','atlocs','niters','smallres',...
        'randseed','voin','vrf_w_std','sim_response_std','parcombs','nm',...
        'npos','ns','-v7.3');
    fprintf('SAVED %s!\n',savefn);
    
    % ~1500 - 2500 seconds per VOI
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function mean_recon = doSimRecon(parcombs,thisdat,nconds,nvox,vrf_w_std,...
    sim_response_std,trn_condList,test_list,trn_all,tst_all,trnX,basis_set,...
    xx,yy,covmat)


% set MATLAB to catch singular matrix warnings as errors
warning('error', 'MATLAB:nearlySingularMatrix');
warning('error', 'MATLAB:singularMatrix');

all_recons = []; all_rpos = []; all_rcond = []; all_parcomb = [];
mean_recon = [];

for p = 1:size(parcombs,1)
    vrf_dat = thisdat;
    
    nochange = find(~parcombs(p,:));
    for ni = 1:length(nochange)
        vrf_dat(:,1:2,nochange(ni)) = repmat(vrf_dat(:,3,nochange(ni)),1,2,1);
    end

    %% generate vRF weights
    for cc = 1:nconds
        for vv = 1:nvox
            % baseline + 2dcos(x,y,sz)*amp
            vrfw(cc,vv,:) = vrf_dat(vv,cc,5) + make2dcos(xx,yy,...
                vrf_dat(vv,cc,1),vrf_dat(vv,cc,2),...
                vrf_dat(vv,cc,3),7)*vrf_dat(vv,cc,4);
        end
    end

    % add noise to the weights
    vrfwn = vrfw + randn(size(vrfw))*vrf_w_std;

    %% generate simulated data from weights

    sim_response = nan(size(trn_condList,1),nvox);
    sim_test = nan(size(test_list,1),nvox);

    % response = sum(channel weight .* channel response to the stimulus)
    % X_all is n_trials x n_channels (with matching channels...)

    for cc = 1:nconds
        tridx = trn_condList==cc;
        sim_response(tridx,:) = trn_all(tridx,:)*squeeze(vrfwn(cc,:,:))';        
        sim_response_n(tridx,:) = sim_response(tridx,:) + mvnrnd(ones(1,nvox)*...
            sim_response_std,covmat{cc},sum(tridx));
        testidx = test_list(:,2)==cc;
        sim_test(testidx,:) = tst_all(testidx,:) * squeeze(vrfwn(cc,:,:))';
        sim_test_n(testidx,:) = sim_test(testidx,:) + mvnrnd(ones(1,nvox)*...
            sim_response_std,covmat{cc},sum(testidx));
    end


    %% compute channel responses

    % make design matrix
    trn = sim_response_n;
    w = trnX\trn;

    for cc = 1:nconds-1
        % apply the EM weight matrix to the test data (L or R)
        tstset = test_list(:,2)==cc;
        tst = sim_test_n(tstset,:);
        try
            x = inv(w*w')*w*tst';
        catch ME
            switch ME.identifier
                case 'MATLAB:nearlySingularMatrix'
                    return;
                case 'MATLAB:singularMatrix'
                    return;
            end
        end
        chan_resp(tstset,:) = x';
        clear x;
    end

    %% weighted sum of basis set = representation

    % basis_set is numChan x res^2
    % chan_resp is nTrials x nChan

    recon_vec = chan_resp*basis_set;

    all_rpos = cat(1,all_rpos,test_list(:,1));
    all_rcond = cat(1,all_rcond,test_list(:,2));
    all_parcomb = cat(1,all_parcomb,repmat(p,size(test_list,1),1));
    all_recons = cat(1,all_recons,recon_vec);

    clear recon_vec w vrfw vrfwn chan_resp x
    clear sim_response sim_response_n sim_test sim_test_n
    clear vrf_dat
end
clear subdat stim tststim trial_list test_list X_all tst_all

%% avg these and calculate error measurements

nstim = length(unique(all_rpos));
clim = [min(all_recons(:)) max(all_recons(:))];

for p = 1:size(parcombs,1)
    for cc = 1:nconds-1
        for ss = 1:nstim
            thisidx = all_rpos==ss & all_rcond==cc & all_parcomb == p;
            mean_recon(p,cc,ss,:) = nanmean(all_recons(thisidx,:));
            clear thisidx;
        end
    end
end