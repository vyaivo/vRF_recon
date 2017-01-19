function simAnalysis_fitSimRecons(niters,voilist,fov,simfix)
%%% This function will fit the reconstructions from the layered IEM
%%% simulations with a 2D cosine. The function inputs are:
%%% 1) niters - number of times to resample the fits
%%% 2) voilist - the regions you are interested in fitting
%%% 3) fov - field of view (in degrees vis. angle) of recon
%%% 4) simfix - a string that modifies the name of the loaded/saved file (to
%%% distinguish different types of layered IEM simulations)
%%%
%%% In order to generate confidence intervals for the figures, the repeated
%%% simulation instances are resampled with replacement before being
%%% averaged and fit.

%% argument parsing
if nargin < 4
    simfix = '_';
%     simfix = '_covNoise_';
end
if nargin < 3
    fov = [13.1250, 9.3750];
end
if nargin < 2
    voilist = {'V1','V2','V3','V3AB','V4','IPS0'};
end
if nargin < 1
    niters = 500;
end

%% some basic params

root = load_root;
simdir = 'sims';

disp('setting up fit params...');
% set up some fit params
s_scales = linspace(0.5,10,20);
FWHM = 1.5625;
all_sp = s_scales .* (FWHM);

fcon = [0,-5; 5,5];
gsteps = [1, 1, all_sp(3)-all_sp(1)];

%% now loop through regions, saving a file each time

skipvoi = cell(length(voilist),1);
for v = 1:length(voilist)
    fnv = sprintf('%s%s/sim%svRFRecon_allSubsNoOutliers_%s.mat',root,...
        simdir,simfix,voilist{v});
    fprintf('Loading %s...\n',fnv);
    load(fnv,'all_sim_recon','skipsub','smallres');
    skipvoi{v} = skipsub;
    fprintf('done.\n');
%     mean_recon = all_sim_recon;
%     clear all_sim_recon;

    if ~exist('nsi','var')
        [nsi,nm,nc,npos,~] = size(all_sim_recon);            
        xx = linspace(-fov(1)/2,fov(1)/2,smallres(1));
        yy = linspace(-fov(2)/2,fov(2)/2,smallres(2));
        [mx, my] = meshgrid(xx, yy);
        evalpts = [reshape(mx,numel(mx),1) reshape(my,numel(my),1)];
    end
    rs2 = 38502462 + v;
    rng(rs2);

    parfor i = 1:niters
        fprintf('iter %d...',i);
        rsi = single(randsample(nsi,nsi,1));
        mr = squeeze(nanmean(all_sim_recon(rsi,:,:,:,:),1));
        all_rs = reshape(mr,nm*nc*npos,[])';

        maxpts = max(all_rs,[],1);
        maxfn = @(x) (find(all_rs(:,x)==maxpts(x)));
        try
            ampis = cell2mat(arrayfun(maxfn, 1:size(all_rs,2),...
                'UniformOutput',0));
        catch
            fprintf('ERROR with MAX...');
            continue;
        end
        [bpar,berr] = gridfit_restrict(all_rs,evalpts,all_sp,...
            evalpts(ampis,:),[0 inf],0,1);

        [bfpar,bferr,~] = gridfit_finetune_con(all_rs,...
            @make2dcos_grid,bpar,evalpts,gsteps,fcon,0,1);

        badfits = find(berr - bferr' < 0);
        for ii = 1:length(badfits)
            recn = badfits(ii);
            % use the original grid fits if finetuning error is high
            bfpar(recn,:) = bpar(recn,:);
            bferr(recn) = berr(recn);
        end

        allfits(i,:,:) = bfpar;
        allerr(i,:) = bferr;
        fprintf('done!\n');
%         clear bpar berr bfpar bferr badfits
    end

    fitfn = sprintf('%s%s/sim%sRMSE_fit_%s.mat',root,simdir,simfix,voilist{v});
    save(fitfn,'allfits','allerr','niters','rs2');
    fprintf('SAVED fits to %s\n',fitfn);
    
    clear allfits allerr
end