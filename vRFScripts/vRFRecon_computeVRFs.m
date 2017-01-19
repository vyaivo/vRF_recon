function vRFRecon_computeVRFs(subj,vois)
% Very simple. Multiplies ridge-regressed weight matrix (from
% vRFRecon_ridgevRFs_wCVThresh.m) by the basis set defined in a pixel grid
% with the FOV/resolution defined below. Saves out a pixel grid image for
% each vRF.

% VAV 3/1/2014
% cleaned up with new file inputs / for OSF. VAV 12/10/16

if nargin < 1
    subj = {'AA','AI','AL','AP','AR','AT','AU'};
end
if nargin < 2
    vois = {'V1','V2','V3','V3AB','V4','IPS0'};
end

root = load_root;

vrf_dir = 'vRFs';
vrf_save_str = 'ridgeVRF_CVBilat';
vrfstr = '_threshVRFRidge_LambdaMinBIC_CVCorr_';
conds = {'Left','Right','Fix'};
nc = length(conds);

% some vRF visualization params (define pixel grid)
fov = [10 5];           % deg, horiz; vert
res = [68 51];

%%
for s = 1:length(subj)
    for v = 1:length(vois)
        
        fn = sprintf('%s%s/%s%s%s.mat',root,vrf_dir,subj{s},...
            vrfstr,vois{v});
        if ~exist(fn,'file')
            fprintf('cannot find file. skipping %s...\n',fn);
            continue;
        else
            fprintf('loading %s...\n',fn);        
            load(fn);
        end
        %%
        if ~exist('mybasis','var')
            
            gridxPts = linspace(-fov(1)/2,fov(1)/2,res(1));
            gridyPts = linspace(-fov(2)/2,fov(2)/2,res(2));
            [gridx, gridy] = meshgrid(gridxPts,gridyPts);

            [cx,cy] = meshgrid(chanX,chanY);
            size_constant = (1/rad2fwhm(1)) * FWHM;
            % build our basis set (want nChannels x nPixels)
            mybasis = build_basis_pts(cx(:),cy(:),size_constant,...
                gridx(:),gridy(:),0);
            mybasis = mybasis.';
        end
                
        vrf_vec = cell(1,nc);
        vrf_mat = cell(1,nc);
        for cc = 1:nc
            % mybasis is nchan x n_recon_pts
            % w_conds is nc x nchan x nvox
            if ismatrix(w_conds)
                % there's only one voxel in this one. no longer need to
                % transpose w_conds
                vrf_vec{cc} = squeeze(w_conds(cc,:))*mybasis;
            else
                vrf_vec{cc} = squeeze(w_conds(cc,:,:))'*mybasis;
            end
            vrf_mat{cc} = nan(res(2),res(1),size(vrf_vec{cc},1));

            % vrf_mat is Y by X by nVoxels
            for t = 1:size(vrf_vec{cc},1)
                vrf_mat{cc}(:,:,t) = reshape(vrf_vec{cc}(t,:),res(2),res(1));
            end
        end
        
        savefn = sprintf('%s%s/%s_%s_%s.mat',root,vrf_dir,subj{s},...
            vrf_save_str,vois{v});
        fprintf('saving to: %s...\n',savefn);
        
        save(savefn,'fn','vrf_vec','vrf_mat','conds','chanX','chanY','FWHM',...
            'size_constant','fov','res','vrfstr','bestvoxhemi');
        
        clear fn savefn
        clear vrf_vec vrf_mat
    end
end