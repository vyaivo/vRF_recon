%%% vRFanalysis_attnModulations_byHemisphere_noOutliers_stats.m
%%% This script does all of the vRF statistical analysis to see how they
%%% change with attention. Here are its component analyses:

%%% 1) Runs a randomized two-way test of attention hemifield (contra vs. ipsi)
%%% by voxel hemipshere (LH vs. RH).

%%% created VAV 6/9/2016; renamed/cleaned for OSF 12/16/2016
%%% VAV 1/2/2017: adapted script to include the hemisphere stuff
%%% VAV 1/13/2017: loaded data from a file of all difference scores,
%%% instead of doing it within this script.
%%% VAV 1/18/2017: made sure all stats were done using subject-level
%%% resampling / bootstrapping

clear all;

%% set up some params

% load the data
root = load_root;
fitdir = 'vRFfits';
fn = sprintf('%s%s/all_vRF_diffscores_noOutliers.mat',root,fitdir);
load(fn);
% fitmat is: data, sub label, VOI label, bin label, distance covariate,
% param label, contra/ipsi (attn side) label, hemisphere label

% output file name
savefn2 = fullfile(root,fitdir,'vRF_attnMods_statsWHemi.mat');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% All analyses WITH the hemisphere/hemifield factor

%% GET MEAN VRF CHANGES
rs2 = sum(clock*10000);
rng(rs2);
niters = 10000;

all_diffscores = nan(npar,nv,na,niters);
mean_diffscores = nan(npar,nv,na);
ci_diffscores = nan(npar,nv,na,2);
pval_diffscores = nan(npar,nv,na);

for p = 1:npar
    for v = 1:nv
        fprintf('Param %d, %s\n',p,voilist{v});
        for a = 1:na
            thisd = fitmat(fitmat(:,3)==v & fitmat(:,6)==p & fitmat(:,7)==a,:);

            parfor i = 1:niters
               % resample across subjects...
               rsi = randsample(ns,ns,1);
               rsd = [];
               for s = 1:ns
                   rsd = [rsd; thisd(thisd(:,2)==rsi(s),1)];
               end
               all_diffscores(p,v,a,i) = mean(rsd);
            end
            tmpscore = squeeze(all_diffscores(p,v,:));
            mean_diffscores(p,v,a) = nanmean(tmpscore);
            ci_diffscores(p,v,a,:) = prctile(tmpscore,[2.5,97.5]);
            pval_diffscores(p,v,a) = 2*min([nanmean(tmpscore<0),nanmean(tmpscore>0)]);

            clear thisd
        end
    end
end
[~,meandiffscore_fdr] = fdr(pval_diffscores,0.05);

%% PLOT MEAN VRF CHANGES
figure;
for p = 1:npar
    subplot(2,2,p);
    ci = squeeze(abs(ci_diffscores(p,:,:,:) - mean_diffscores(p,:,:)));
    [h1, h2] = barwitherr(ci,squeeze(mean_diffscores(p,:,:)));
    h1(1).Parent.XTickLabel = voilist;
    h1(1).Parent.XTickLabelRotation = -35;
%     
%     sigv = find(hemifdr(p,:));
%     if ~isempty(sigv)
%         hold on;
%         plot(sigv,repmat(h1(1).Parent.YLim(2)*0.9,1,length(sigv)),'k*');
%     end
    
    box off;
    xlim([0 nv+1]);
    title(parlabel{p});
end
legend('contra','ipsi');

%% 'bootstrapped' ANOVA
% e.g. just do main effects (mean A1 - mean A2) <> 0, etc.
% and interaction (mean A1B1 - mean A1B2) - (mean A2B1 - mean A2B2) <> 0

rs4 = 80459740;
rng(rs4);

for p = 1:npar
    fprintf('Param %d\n',p);
    for v = 1:nv
        thisdat = fitmat(fitmat(:,3)==v & fitmat(:,6)==p,:);
        
        iter_attn = nan(niters,1);
        iter_hemi = nan(niters,1);
        iter_attnbyhemi = nan(niters,1);
        
        parfor i = 1:niters
           % resample across subjects...
           rsi = randsample(ns,ns,1);
           rsd = [];
           for s = 1:ns
               rsd = [rsd; thisdat(thisdat(:,2)==rsi(s),:)];
           end
            % main effect of attn hemifield
            iter_attn(i) = nanmean(rsd(rsd(:,7)==1,1)) - ...
                nanmean(rsd(rsd(:,7)==2,1));

            % main effect of vox hemisphere
            iter_hemi(i) = nanmean(rsd(rsd(:,8)==1,1)) - ...
                nanmean(rsd(rsd(:,8)==2,1));
    
            % interaction effect
            % (mean A1B1 - mean A1B2) - (mean A2B1 - mean A2B2) <> 0
            iter_attnbyhemi(i) = ( nanmean(rsd(rsd(:,7)==1 & rsd(:,8)==1,1)) - ...
                nanmean(rsd(rsd(:,7)==1 & rsd(:,8)==2,1)) ) - ...
                ( nanmean(rsd(rsd(:,7)==2 & rsd(:,8)==1,1)) - ...
                nanmean(rsd(rsd(:,7)==2 & rsd(:,8)==2,1)) );
        end
        
        p_attn(p,v) = 2*min( [nanmean(iter_attn<0),nanmean(iter_attn>0)] );
        p_hemi(p,v) = 2*min( [nanmean(iter_hemi<0),nanmean(iter_hemi>0)] );
        p_attnbyhemi(p,v) = 2*min( [nanmean(iter_attnbyhemi<0),...
            nanmean(iter_attnbyhemi>0)] );

        clear thisdat
    end
    
    [pfdr_attn(p), pmask_attn(p,:)] = fdr(p_attn(p,:), 0.05);
    [pfdr_hemi(p), pmask_hemi(p,:)] = fdr(p_hemi(p,:), 0.05);
    [pfdr_attnbyhemi(p), pmask_attnbyhemi(p,:)] = fdr(p_attnbyhemi(p,:), 0.05);
end

%% BOOTSTRAP CIs ON ALL GROUPS
rs5 = 490249063;
rng(rs5);

for p = 1:npar
    for v = 1:nv
        thisd = fitmat(fitmat(:,3)==v & fitmat(:,6)==p,:);
        
        cl = nan(1,niters); cr = cl; il = cl; ir = cl;
        parfor i = 1:niters
            % resample across subjects...
            rsi = randsample(ns,ns,1);
            rsd = [];
            for s = 1:ns
                rsd = cat(1,rsd,thisd(thisd(:,2)==rsi(s),:));
            end
            cl(i) = nanmean(rsd(rsd(:,7)==1&rsd(:,8)==1,1));      
            cr(i) = nanmean(rsd(rsd(:,7)==1&rsd(:,8)==2,1));  
            il(i) = nanmean(rsd(rsd(:,7)==2&rsd(:,8)==1,1));  
            ir(i) = nanmean(rsd(rsd(:,7)==2&rsd(:,8)==2,1));
        end
        all_groupbootmean(p,v,:) = cellfun(@nanmean,{cl,cr,il,ir});
        all_groupbootci(p,v,1,:) = abs(cellfun(@prctile,{cl,cr,il,ir},...
            repmat({2.5},1,4)) - squeeze(all_groupbootmean(p,v,:))' );
        all_groupbootci(p,v,2,:) = abs(cellfun(@prctile,{cl,cr,il,ir},...
            repmat({97.5},1,4)) - squeeze(all_groupbootmean(p,v,:))' );
    end
end

%% PLOT GROUP MEANS AS LINES

% yl = [-1 1; -1 1; 0 0.5; -0.2 0.2];
yl = [-2 0.5; -1 0.75; 0 1; -0.2 0.2];

figure; pp=1;
for p = 1:npar
    for v = 1:nv
        groupci = squeeze(all_groupbootci(p,v,:,:))';
        groupmean = squeeze(all_groupbootmean(p,v,:));
%         
        hh(v) = subplot(npar,nv,pp); hold all;
        h = errorbar(1:2,groupmean(1:2),groupci(1:2,1),groupci(1:2,2),...
            'LineWidth',1.5);
        errorbar(1:2,groupmean(3:4),groupci(3:4,1),groupci(3:4,2),...
            'LineWidth',1.5);

%         h = plot(1:2,groupmean(1:2),'-o','LineWidth',2,'MarkerSize',5);
%         plot(1:2,groupmean(3:4),'-o','LineWidth',2,'MarkerSize',5);
        
        if pmask_attn(p,v) == 1
            % * is main effect of attn
            plot(1.1, 0.9*yl(p,2),'k*');
        end
        if pmask_hemi(p,v) == 1
            % o is main effect of hemisphere
            plot(1.5, 0.9*yl(p,2),'ko');
        end
        if pmask_attnbyhemi(p,v) == 1
            % x is interaction effect
            plot(1.9, 0.9*yl(p,2),'kx');
        end
        
        h.Parent.XTick = 1:2;
        if p == 4
            h.Parent.XTickLabel = {'LH','RH'};
            h.Parent.XTickLabelRotation = -40;
        else
            h.Parent.XTickLabel = '';
        end
        if pp == 1
            legend('contra','ipsi','attn hemifield','vox hemisphere','interaction');
        end
        box off;
        xlim([0 3]);
        title([parlabel{p},', ',voilist{v}]);
        
        pp = pp+1;
        
    end
    match_ylim(hh,yl(p,:));
    
    clear hh
end


%% SAVE

save(savefn2,'rs2','all_diffscores','mean_diffscores',...
    'ci_diffscores','pval_diffscores','allout','rs3','pval_ci2',...
    'hemifdr','rs4','p_attn','p_hemi',...
    'p_attnbyhemi','pfdr_attn','pmask_attn','pfdr_hemi','pmask_hemi',...
    'pfdr_attnbyhemi','pmask_attnbyhemi','-v7.3');
fprintf('SAVED %s\n',savefn2);
