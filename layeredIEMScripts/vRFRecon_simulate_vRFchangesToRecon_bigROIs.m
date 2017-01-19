function vRFRecon_simulate_vRFchangesToRecon_bigROIs(doROI)
% giant script to run through all the possible / desired combos of vrf
% changes. uses compiled vRF data to exclude simulations of voxels with
% difference scores > 3*SD of the population mean.

% VAV 12/22/2016
% VAV 1/6/2017 -- added this for big ROIs

%%
if nargin < 1
    doROI = 1;
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
vrfn = fullfile(root,fdir,'AllSubs_CompiledvRFs_bigROIs.mat');
load(vrfn);

rfn = sprintf('%s%s/AllSubs_CompiledRecons_bigROIs_%s.mat',root,'reconfits',...
    'goodvRFs_TrnAllJitter');
rd = load(rfn,'skipsub');

% % get vRF outliers (based on vRF diff scores > 3*sd from mean)
% voutfn = fullfile(root,fdir,'vRF_attnMods_stats.mat');
% load(voutfn,'allout');

% define some params based on the first sub/VOI
exampleFile = sprintf('%s%s/%s_goodvRFs_TrnAllJitter_acrossSess_superVis.mat',...
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
    
    savefn = sprintf('%s%s/sim_vRFRecon_allSubsNoOutliers_%s.mat',root,simdir,...
        voilist{voin});
    tic
    for sn = 1:ns
        %%
        sub = sublist{sn};
        subdat = alldat{sn,voin};
        
        % need to convert vRF FWHM size to size constant for make2dcos fcn
        subdat(:,:,3) = subdat(:,:,3) ./ rad2fwhm(1);

        if isempty(subdat) || size(subdat,1) == 1 || ismember(sn,skipsub) || ...
                ismember(sn,rd.skipsub{voin})
            % if there's not enough data or we should skip the subject, skip
            if ~ismember(sn,skipsub)
                skipsub = [skipsub sn];
            end
            fprintf('SKIPPING %s %s, not enough vox\n',sub,voilist{voin});
            continue;
        end
        nvox = size(subdat,1);

        %%
        
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
            disp(ii);
            tmp = doSimRecon(parcombs,subdat,nconds,...
                nvox,vrf_w_std,...
                sim_response_std,trn_condList,test_list,trn_all,tst_all,...
                trnX,basis_set,xx,yy);
            if isempty(tmp)
                skipsub = [skipsub,sn];
                continue;
            end
            iterrecon(ii,:,:,:,:) = tmp;
        end
        subrecon(sn,:,:,:,:,:) = iterrecon;
    end
    toc
    all_sim_recon = squeeze(nanmean(subrecon,1));
    save(savefn,'all_sim_recon','skipsub','atlocs','niters','smallres',...
        'randseed','voin','vrf_w_std','sim_response_std','parcombs','nm',...
        'npos','ns','-v7.3');
    fprintf('SAVED %s!\n',savefn);
    
    % ~2500 seconds per VOI
end
