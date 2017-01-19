% compile_vrf_data_allsubs.m
% vav 5/28/2015
% just creates a big matrix file of all vrf params, sub & VOI labels
% this can be used for the simulation analysis. doesn't contain bin labels
% cleaned up for OSF 12/12/2016

atlocs = [-2.1074, 0.3664; 2.1233, 0.3664];
plotbool = 1;
sublist = {'AA','AI','AL','AP','AR','AT','AU'};
voilist = {'V1','V2','V3','V3AB','V4','IPS0'};

%% load the data
root = load_root;
fitdir = 'vRFfits';
rfdir = 'vRFs';

npar = 5;               % x,y,size,amp,baseline
nc = 3;
ns = length(sublist);
nv = length(voilist);
alldat = cell(ns,nv);
all_rmse = [];
saveWHemi = 1;

%%
for v = 1:length(voilist)
    
    threshold_stats(v).nvox = 0;
    threshold_stats(v).thresh1 = 0;
    threshold_stats(v).thresh2 = 0;
    threshold_stats(v).thresh3 = 0;    
    
    voi_rmse = [];
    
    for s = 1:length(sublist)
        %% load fit stuff
        fn = sprintf(...
            '%s%s/%s_%s_GridFit_ThreshVRFs_LambdaMinBIC_ALLSESS.mat',...
            root,fitdir,sublist{s},voilist{v});
        if ~exist(fn,'file')
            disp(['Skipping ' sublist{s} ', ' voilist{v}]);
            continue;
        end
        load(fn,'bfpar','bferr');
        voi_rmse = cat(1,voi_rmse,bferr{3});
        
        %% load the vox indices that we want to keep
        fn2 = sprintf(...
            '%svRFs/%s_threshVRFRidge_LambdaMinBIC_CVCorr_%s.mat',...
            root,sublist{s},voilist{v});
        load(fn2,'thresh_r2_vox','thresh_ridge_vox','thresh_cv_vox',...
            'ntotalvox','bestvox','bestvoxhemi');
        threshold_stats(v).nvox = threshold_stats(v).nvox + ntotalvox;
        threshold_stats(v).thresh1 = threshold_stats(v).thresh1 + numel(thresh_r2_vox);
        threshold_stats(v).thresh2 = threshold_stats(v).thresh2 + numel(thresh_ridge_vox);
        threshold_stats(v).thresh3 = threshold_stats(v).thresh3 + numel(thresh_cv_vox);
        
        %% put everything in an array
        for c = 1:nc
            newdat(:,c,:) = bfpar{c};
            suberr(:,c,:) = bferr{c};
        end
        subdat = cat(3,newdat,nan(size(newdat,1),nc,2));
        all_hemilabs{s,v} = bestvoxhemi;        
        
        % also calc fix dist from attend location        
        xc = subdat(:,3,1);
        yc = subdat(:,3,2);
        for c = 1:nc-1
            % rho (for polar coord bins)
            subdat(:,c,npar+1) = sqrt( (xc-atlocs(c,1)).^2 + (yc-atlocs(c,2)).^2 );
            % theta (for polar coord bins)
            subdat(:,c,npar+2) = atan2( yc-atlocs(c,2), xc-atlocs(c,1) );
        end

        % convert size constant to FWHM
        subdat(:,:,3) = rad2fwhm(subdat(:,:,3));
        
        % save out all data in cell structs
        alldat{s,v} = subdat;
        allerr{s,v} = suberr;
        allgoodvox{s,v} = bestvox;

        clearvars newdat subdat suberr bfpar bferr bestvox
    end
    
    threshold_stats(v).percThresh = threshold_stats(v).thresh3 ./ threshold_stats(v).nvox;
    threshold_stats(v).rmse = mean(voi_rmse);
    all_rmse = cat(1,all_rmse,voi_rmse);
end

vrf_threshhold_table = struct2table(threshold_stats,'RowNames',voilist);

%% save out this stuff for later use

savefn = fullfile(root,fitdir,'AllSubs_CompiledvRFs.mat');
save(savefn,'alldat','allerr','allgoodvox','sublist','voilist','atlocs',...
    'npar','nc','ns','nv','vrf_threshhold_table','all_hemilabs');
fprintf('Saved %s\n',savefn);