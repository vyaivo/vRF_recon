function vRFRecon_gridfit_maxRestrict_subRecons(subname,voi_names,trnPrefix)
% Fits per subject reconstructions for a set of VOIs. The center of the fit
% is restricted to an area near the maximum of the recon (right now, +/- 1
% degrees of visual angle).
% utilizes TCS gridfit toolbox (gridfit.m)

% VAV 11/1/14
% last edited 9/29/2016

root = load_root;
fitdir = 'reconfits';

if nargin < 3
    trnPrefix = 'TrnAllJitter';
end
analysisPrefix = 'acrossSess';

if ~exist(fullfile(root,fitdir),'dir')
    mkdir(fullfile(root,fitdir));
end

% The average recon file contains the following variables: subs, voi_names,
% voi_avg (cell array of VOIs containing cond x position x Xdim x Ydim)
% and, if enabled, voiflipavg (cell array of VOIs containing position x
% Xdim x Ydim).

for v = 1:length(voi_names)
    
    fname = sprintf('%srecons/%s_%s_%s_%s.mat',...
        root,subname,trnPrefix,analysisPrefix,voi_names{v});    
    if ~exist(fname)
        disp(['Cannot find ' fname ', skipping...']);
        return;
    end
    load(fname);
    %% Make parameter grid
    if v == 1      
        xx = linspace(-fov(1)/2,fov(1)/2,res(1));
        yy = linspace(-fov(2)/2,fov(2)/2,res(2));
        [mx, my] = meshgrid(xx,yy); clear xx; clear yy;
        evalpts = [reshape(mx,numel(mx),1) reshape(my,numel(my),1)];
        % grid parameters        
        s_scales = linspace(0.5,10,20);
        allp = s_scales .* (FWHM);	% this is now the scaled size constant (s)
                                    % to convert to dva FWHM, do rad2fwhm
    end
    
    fns = sprintf('%s%s/%s_%s_ReconGridFits_%s_%s.mat',...
        root,fitdir,subname,voi_names{v},analysisPrefix,trnPrefix);
    disp(['Fitting ' voi_names{v}]);
    [bfpar, bferr] = doFit(recon_avg,allp,evalpts,1);
    
    save(fns,'bfpar','bferr','recon_avg','evalpts','res','locs','-v7.3');
    disp(['Saved ' fns]);
    
    clear bfpar bferr recon_avg
end

function [bfpar, bferr] = doFit(recon_vec,allp,evalpts,parflag)
% give fmincon some parameter constraints -- e.g., x,y not out of
% bounds; size not bigger than screen [LB;UB]
fcon = [0,-5; 5,5];
gsteps = [1, 1, allp(2)-allp(1)];

for cc = 1:size(recon_vec,1)
    %% make grids
    disp(['Fitting condition ' num2str(cc)]);
    data = squeeze(recon_vec(cc,:,:))';
    [rinds,~] = find(data==repmat(max(data),size(data,1),1));
    xp = evalpts(rinds,1);
    yp = evalpts(rinds,2);
    [bpar,berr] = gridfit_restrict(data,evalpts,allp,[xp yp],[0 inf],parflag);
    
    %% grid fit & fmincon fine tuning
    disp(['Doing recon grid fit for condition ' num2str(cc),'...']);
    [bfpar(cc,:,:), bferr(cc,:), ~] = gridfit_finetune_con(data,...
        @make2dcos_grid,bpar,evalpts,gsteps,fcon,parflag);
    badfits = find(berr - squeeze(bferr(cc,:)) < 0);
    for i = 1:length(badfits)
        recn = badfits(i);
        % use the original grid fits, for now
        bfpar(cc,recn,:) = bpar(recn,:);
        bferr(cc,recn) = berr(recn);
    end
    clear badfits;
end