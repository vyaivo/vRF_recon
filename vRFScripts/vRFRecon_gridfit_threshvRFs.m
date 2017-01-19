function vRFRecon_gridfit_threshvRFs(subj,voi_names)
% utilizes a gridfit toolbox -- make sure everything is in path (gridfit.m)
% VAV 10/26/14
% cleaned up file/var names for OSF. VAV 12/12/16

if nargin < 1
    subj = {'AA','AI','AL','AP','AR','AT','AU'};
end
if nargin < 2
    voi_names = {'V1','V2','V3','V3AB','V4','IPS0'};
end
%%
root = load_root;
vrfdir = 'vRFs';
vrfn = 'ridgeVRF_CVBilat';
fitdir = 'vRFfits';

if ~exist(fullfile(root,fitdir),'dir')
    mkdir(fullfile(root,fitdir));
end

% gridfit params
xy_density = 1; % number of samples per dva for mux/muy spacing
szlist = logspace(0.01,1,20);    % FWHM for search grid

% if parpool not already started, start one....
p = gcp('nocreate');
if isempty(p)
    parpool(12);
end

%%
for s = 1:length(subj)
    for v = 1:length(voi_names)
        fprintf('Fitting %s, %s \n',subj{s},voi_names{v});        
        fn = sprintf('%s%s/%s_%s_%s.mat',root,vrfdir,subj{s},vrfn,voi_names{v});
        if ~exist(fn,'file')
            fprintf('No data found for %s, %s. Skipping...\n', subj{s}, ...
                voi_names{v});
            continue;
        else            
            fprintf('loading %s...\n',fn);
        end
        load(fn);
        
        % figure out where to evaluate the grid
        xx = linspace(-fov(1)/2,fov(1)/2,res(1));
        yy = linspace(-fov(2)/2,fov(2)/2,res(2));
        [mx, my] = meshgrid(xx,yy); clear xx; clear yy;
        evalpts = [reshape(mx,numel(mx),1) reshape(my,numel(my),1)];
        % grid parameters
        xgrid = linspace((-fov(1)/2)+0.5,(fov(1)/2)-0.5,xy_density*fov(1)+1);
        ygrid = linspace((-fov(2)/2)+0.5,(fov(2)/2)-0.5,xy_density*fov(2)+1);
        szgrid = szlist ./ rad2fwhm(1);
        [gp1, gp2, gp3] = ndgrid(xgrid,ygrid,szgrid);
        gsteps = [xgrid(3)-xgrid(1), ygrid(3)-ygrid(1), szgrid(3)-szgrid(1)];

        disp('Making search grid...')
        grid_params = [gp1(:) gp2(:) gp3(:)];
        
        %% make the actual grid for input to the function
        mygrid = make_grid(@make2dcos_grid,evalpts,grid_params);
        fcon = [0,-5; 5,5];
        for cc = 1:length(vrf_vec)
            % gridfit function needs to have voxels in the 2nd dimension
            data = vrf_vec{cc}';
            
            fprintf('Running grid fit for condition: %s...\n', conds{cc});
            [bpar,berr,bfcn] = gridfit(data,mygrid,grid_params,[0 inf]);
            [bpar2, berr2, bfcn2] = gridfit_finetune_con(data,@make2dcos_grid,...
                bpar,evalpts,gsteps,fcon);
            % if the finetuning doesn't result in a better fit, use the
            % best grid fit.
            worsefit = find(berr' < berr2);
            bpar2(worsefit) = bpar(worsefit);
            berr2(worsefit) = berr(worsefit);
            bfcn2(worsefit) = bfcn(worsefit);
            bfpar{cc} = bpar2;
            bferr{cc} = berr2;
            bffcn{cc} = bfcn2;
            clear bpar berr bfcn bpar2 berr2 bfcn2
        end
        fns = sprintf('%s%s/%s_%s_GridFit_ThreshVRFs_LambdaMinBIC_ALLSESS.mat',...
            root,fitdir,subj{s},voi_names{v});
        fprintf('Saving grid fits to %s\n', fns);
        save(fns,'bfpar','bferr','bffcn','vrf_vec','worsefit','res',...
            'evalpts','bestvoxhemi','-v7.3');
    end
end

% delete(p);