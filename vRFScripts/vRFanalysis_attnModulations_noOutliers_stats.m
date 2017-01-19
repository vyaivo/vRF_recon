%%% vRFanalysis_attnModulations_noOutliers_stats.m
%%% This script does all of the vRF statistical analysis to see how they
%%% change with attention. Here are its component analyses:

%%% 1) Calculates difference scores for each vRF fit parameter (e.g.,
%%% attend - fixation values).
%%% 2) Finds difference score outliers (3*SD > mean).
%%% 3) Bootstraps the mean vRF difference scores in each ROI
%%% 4) Fits polynomials of order 0, 1, & 2 (e.g. flat line, sloped line,
%%% and quadratic) to the vRF difference scores as a function of the
%%% voxel's distance from the attention target during the attend fixation
%%% condition. Cross-validates these fits to find the best ones.
%%% VAV 1/2/2017: this script was adapted to attach a contralateral vs
%%% ipsilateral label to each voxel, and to repeat analyses 3 & 4 with an
%%% additional factor.

%%% VAV 6/9/2016; edited for OSF 12/16/2016
clear all;

%% set up some params

% load the data
root = load_root;
fitdir = 'vRFfits';
% fn = fullfile(root,fitdir,'AllSubs_CompiledvRFs.mat');
fn = sprintf('%s%s/all_vRF_diffscores_noOutliers.mat',root,fitdir);
load(fn);

parlabel = {'position','size','amplitude','baseline'};
% analysisType will set a switch to do the vRF analyses with or without the
% attention hemifield factor (e.g., contra or ipsi)

% save file names are different for each analysis type
savefn1 = fullfile(root,fitdir,'vRF_attnMods_stats.mat');

% %% sort the vRFs
% 
% % This makes a giant matrix with vRF difference scores as the first column,
% % and all subsequent columns are labels that describe properties of the
% % vRF. i.e.
% 
% % col 1: data
% % col 2: subject
% % col 3: region of interest (or VOI, volume of interest)
% % col 4: bin number (by vRF's distance from the attended target)
% % col 5: distance from the attended target (a continuous value)
% % col 6: parameter label (position, size, amplitude, baseline)
% % col 7: contralateral vs. ipsilateral to attention target
% 
% % set up dist from attend bins
% distbins = 0:0.25:2.5;
% 
% % create empty label vectors to fill in
% adat = [];slab = [];vlab = [];dlab = [];
% blab = [];clab = [];voxlab =[];hlab=[];
% 
% for s = 1:ns
%     for v = 1:nv
%         % read in data from each sub/VOI
%         if isempty(alldat{s,v})
%             continue;
%         end
%         %%
%         thisdat = alldat{s,v}(:,:,1:npar);
%         
%         %% calculate vRF param changes (attend - fix difference scores)
%         
%         attndat = nan(size(thisdat,1),2,size(thisdat,3)-1);
%         contra_ipsi = nan(2,size(thisdat,1));
%         for c = 1:nc-1
%             
%             % distance from attention locus during fix cond
%             dist_at_fix(c,:) = sqrt( (thisdat(:,3,1)-atlocs(c,1)).^2 + ...
%                 (thisdat(:,3,2)-atlocs(c,2)).^2 );            
%             % sort vRFs by distance from attention locus during fix cond
%             [~,bi(c,:)] = histc(dist_at_fix(c,:),distbins);
%             
%             % distance from attention locus during attention cond
%             dist_attend = sqrt( (thisdat(:,c,1)-atlocs(c,1)).^2 + ...
%                 (thisdat(:,c,2)-atlocs(c,2)).^2 );
%             
%             % now calc diff between position measures above
%             attndat(:,c,1) = dist_attend' - dist_at_fix(c,:);
%             
%             for p = 3:npar
%                 % calculate difference from fixation condition
%                 attndat(:,c,p-1) = thisdat(:,c,p) - thisdat(:,3,p);
%             end
%             
%             % LABEL CONTRA VS IPSI
%             % A) first get hemi labels
%             % B) then get attn labels. if they ARE OPPOSITE (e.g., hemi = 1
%             % & attn = 2), you get CONTRA (1). if they ARE THE SAME you get
%             % IPSI (2).
%             ipsi = find(all_hemilabs{s,v} == c);
%             contra = setxor(ipsi,1:size(thisdat,1));
%             
%             contra_ipsi(c,contra) = 1;
%             contra_ipsi(c,ipsi) = 2;
%             
%             clear dist_attend
%         end
%         
%         vdat = repmat([1:size(attndat,1)]',1,2,npar-1);
%         hemidat = repmat(all_hemilabs{s,v},1,2,npar-1);
%         %%
%         % we want to fold the data together across conditions (attend L vs R)
%         % so we should get 2x the vRFs
%         catdat = squeeze(cat(1,attndat(:,1,:),attndat(:,2,:)));
%         cathemi = squeeze(cat(1,hemidat(:,1,:),hemidat(:,2,:)));
%         catvoxdat = squeeze(cat(1,vdat(:,1,:),vdat(:,2,:)));
%         catbin = [bi(1,:) bi(2,:)];
%         catfdist = [dist_at_fix(1,:), dist_at_fix(2,:)];
%         catcipsi = [contra_ipsi(1,:),contra_ipsi(2,:)];
%         
%         % thisdat has size nvRF x npar
%         sl = repmat(s,size(catdat));                % subject label
%         vl = repmat(v,size(catdat));                % VOI label
%         cl = repmat(catcipsi',1,size(catdat,2));     % contra/ipsi label
%         bl = repmat(catbin',1,size(catdat,2));      % bin label
%         dl = repmat(catfdist',1,size(catdat,2));    % distance fr attend
%         
%         % save out all the labels for the ANOVA
%         slab = cat(1,slab,sl);
%         vlab = cat(1,vlab,vl);
%         blab = cat(1,blab,bl);
%         dlab = cat(1,dlab,dl);
%         clab = cat(1,clab,cl);
%         hlab = cat(1,hlab,cathemi);
%         adat = cat(1,adat,catdat);
%         voxlab = cat(1,voxlab,catvoxdat);
%                 
%         clearvars dist_at_fix catdat catbin catfdist catcipsi vdat catvoxdat
%         clearvars newdat attndat bi sl vl bl dl cl
%     end
% end
% 
% % update number of params, now that we've collapsed the x/y fit parameters
% % into an overarching 'position' (or distance fr. attn targ) metric
% npar = npar-1;
% % also set the number of attn hemifield labels (contra vs ipsi)
% na = numel(unique(clab(:)));
% 
% % set up parameter labels
% plab = nan(size(adat));
% for p = 1:npar
%     plab(:,p) = p;
% end
% 
% % create the giant data matrix!
% fitmat = [adat(:) slab(:) vlab(:) blab(:) dlab(:) plab(:) clab(:) hlab(:)];
% % fitmat is: data, sub label, VOI label, bin label, distance covariate,
% % param label, contra/ipsi (attn side) label, hemisphere label
% 
% clearvars adat slab vlab blab dlab plab clab hlab
% clearvars alldat all_hemilabs
% 
% voxlab2 = voxlab(:);
% 
% %% Get outlier indices if vox changes > 3SD from mean
% % Since the outlier thresholding is based on the difference scores across
% % the whole population (e.g., all conditions, all subjects, both
% % hemispheres), we need to first find the outliers, then save out
% % identifiers which specify which subject/region of interest those outliers
% % belong to. Those indices are saved out in the final analysis data file so
% % they can be used to constrain subsequent analyses.
% 
% outi = cell(nv,1);
% allout = cell(ns,nv);
% for v = 1:nv
%     outi{v} = [];
%     for p = 1:npar
%         % vRF diff scores should be similar to each other in each fit
%         % parameter (e.g., size or amplitude) and ROI. So grab all the data
%         % that fits those categories...
%         rowi = find(fitmat(:,3)==v & fitmat(:,6)==p);   % row indices from fitmat
%         thisdat = fitmat(rowi,1);                       % actual diff scores
%         % find standard deviation of this population
%         sd = std(thisdat);
%         sout = find(abs(thisdat - mean(thisdat)) > 3*sd);
%         % save fitmat indices to outliers
%         outi{v} = cat(1,outi{v},rowi(sout));
%         
%         clear thisdat sd sout rowi
%     end
%     
%     % want to save out voxel outlier indices to allout{s,v}. 
%     % i.e., figure out the subject labels for each of the outliers in voilabel == v
%     
%     % find unique sub labels for these voxels
%     thissubs = unique(fitmat(outi{v},2));
%     % now for each subject, need to get the original voxel # label (e.g.,
%     % vox 35 for sub 1, VOI 1; vox 9 for sub 2, VOI 1, etc.)
%     for si = 1:length(thissubs)
%         % grab fitmat row indices for this sub/VOI
%         thiss = find(fitmat(:,2)==thissubs(si) & fitmat(:,3)==v);
%         % now see if any of them match the voxel outliers we found (i.e. a
%         % set intersection of all voxels for this sub/VOI and for all
%         % outlier voxels for this VOI)
%         subouti = intersect(outi{v},thiss);
%         % now get the specific voxel #s for this sub/VOI and save them
%         allout{thissubs(si),v} = voxlab2(subouti);
%         
%         clear thiss subouti
%     end
%     
% end
% 
% %% remove the outliers for these analyses
% newfmat = fitmat;
% allouti = vertcat(outi{:});     % this is all the outliers in the giant matrix fitmat
% newfmat(allouti,:) = [];        % take them out!
% 
% fitmat = newfmat;
% clear newfmat allouti;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% All analyses without the attn hemifield factor

%% MEAN vRF CHANGES

rs2 = sum(clock*10000);
rng(rs2);
niters = 10000;

all_diffscores = nan(npar,nv,niters);
mean_diffscores = nan(npar,nv);
ci_diffscores = nan(npar,nv,2);
pval_diffscores = nan(npar,nv);

fprintf('CALCULATING MEAN vRF CHANGES...');

for p = 1:npar
    for v = 1:nv
        thisd = fitmat(fitmat(:,3)==v & fitmat(:,6)==p,:);

        parfor i = 1:niters
           % resample across subjects...
           rsi = randsample(ns,ns,1);
           rsd = [];
           for s = 1:ns
               rsd = [rsd; thisd(thisd(:,2)==rsi(s),1)];
           end
           all_diffscores(p,v,i) = mean(rsd);
        end
        tmpscore = squeeze(all_diffscores(p,v,:));
        mean_diffscores(p,v) = nanmean(tmpscore);
        ci_diffscores(p,v,:) = prctile(tmpscore,[2.5,97.5]);
        pval_diffscores(p,v) = 2*min([nanmean(tmpscore<0),nanmean(tmpscore>0)]);

        clear thisd
    end
end

[~,meandiffscore_fdr] = fdr(pval_diffscores,0.05);

%% PLOT MEAN VRF CHANGES

figure;
for p = 1:npar
    subplot(2,2,p);
    ci = abs(ci_diffscores(p,:,:) - mean_diffscores(p,:));
    h(p) = barwitherr(ci,mean_diffscores(p,:));
    h(p).Parent.XTickLabel = voilist;
    h(p).Parent.XTickLabelRotation = -35;
    box off;
    xlim([0 nv+1]);
    title(parlabel{p});
end
drawnow;
fprintf('done\n');

%% FIT PARAM RELATIONSHIPS
rs = sum(clock*10000);
rng(rs);
niters = 10000;

coefs = cell(npar,nv,3);
allssres = nan(npar,nv,niters,3);
allr2s = nan(npar,nv,niters,3);

bestmods = nan(npar,nv);
bestr2s = nan(npar,nv);
allFvals = nan(npar,nv,2);
allpvals = nan(npar,nv,2);

% set this up to catch errors
warning('error', 'MATLAB:rankDeficientMatrix');
warning('error', 'MATLAB:polyfit:PolyNotUnique');

fprintf('FITTING POLYNOMIALS TO vRF DIFF SCORES...');
for p = 1:npar
    for v = 1:nv

        %% first fit each model using 50/50 cross-validation
        thisd = fitmat(fitmat(:,3)==v & fitmat(:,6)==p,:);

        db = thisd(:,5);
        dat = thisd(:,1);
%                 nanb = isnan(dat);
%                 dat(nanb) = [];
%                 db(nanb) = [];

        ssres = nan(niters,3); sstot = nan(niters,3); r2 = nan(niters,3);
        trncv = logical(zeros(numel(dat),niters)); tstcv = trncv;
        parfor i = 1:niters
            [trncv(:,i),tstcv(:,i)] = crossvalind('HoldOut',numel(dat),0.3);
        end

        for n = 0:2
            sr = nan(niters,1); st = sr; rr = sr;
            allf = nan(niters,n+1);
            parfor ii = 1:niters
                trn = trncv(:,ii);
                tst = tstcv(:,ii);
                fpar = polyfit(db(trn),dat(trn),n);
                feval = polyval(fpar,db(tst));
                sr(ii) = sum( (dat(tst) - feval).^2 );
                st(ii) = sum( (dat(tst) - mean(dat(tst))).^2 );
                rr(ii) = 1 - (sr(ii)/st(ii));
                allf(ii,:) = fpar;
            end
            ssres(:,n+1) = sr;
            sstot(:,n+1) = st;
            r2(:,n+1) = rr;
            coefs{p,v,n+1} = allf;
        end

        allssres(p,v,:,:) = ssres;
        allr2s(p,v,:,:) = r2;

        total_ssres = mean(ssres,1);

        % now do a nested F-tests to determine the best model
        for nm = 1:2
            F(nm) = ( (total_ssres(nm)-total_ssres(nm+1)) / 1 ) / ...
                    ( total_ssres(nm+1)/(numel(dat)-(nm+1)) );
            pp(nm) = 1 - fcdf(F(nm),1,numel(dat)-2);
        end
        allFvals(p,v,:) = F;
        allpvals(p,v,:) = pp;

        clear F pp r2
    end
end

% get p-values for best fit coeffs
coef_pval = nan(npar,nv);
bestcoefs_mean = cell(npar,nv);
bestcoefs_ci = cell(npar,nv);

for p = 1:npar
    for v = 1:nv
        % These are comparing n = 1 > n = 0 & n = 2 > n = 1. So if
        % neither of these is the best, default to n = 0 (e.g.,
        % bestmod = 1).
        bestmod = max(find(squeeze(allpvals(p,v,:)) < .05)) + 1;
        if isempty(bestmod)
            bestmod = 1;
        end
        bestmods(p,v) = bestmod;
        bestr2s(p,v) = mean(allr2s(p,v,:,bestmod),3);
        bestcoefs_mean{p,v} = mean(coefs{p,v,bestmod},1);
        bestcoefs_ci{p,v} = prctile(coefs{p,v,bestmod},[2.5, 97.5]);
        coef_pval(p,v) = 2*min([mean(coefs{p,v,bestmod}<0),...
            mean(coefs{p,v,bestmod}>0)]);

    end
end
[fdrthresh,fdrmask] = fdr(coef_pval,0.05);
fprintf('done. Saving... \n');

tab2 = cell2table(bestcoefs_mean','VariableNames',...
    {'Position','Size','Amplitude','Baseline'},'RowNames',voilist);

%% save
save(savefn1,'rs','allssres','allr2s','bestmods','bestr2s','bestcoefs_mean',...
    'bestcoefs_ci','allFvals','allpvals','niters','coefs','coef_pval',...
    'fdrthresh','fdrmask','rs2','all_diffscores','mean_diffscores',...
    'ci_diffscores','pval_diffscores','allout','tab2','-v7.3');
fprintf('SAVED %s\n',savefn1);

%% plot best poly fits for all ROI/params        
cm = lines(nv);
figure; pp = 1;
for p = 1:npar
    for v = 1:nv
        h(pp) = subplot(npar,nv,pp);

        % plot scatter (each dot = 1 vRF)
        thisx = fitmat(fitmat(:,3)==v & fitmat(:,6)==p,5);
        h2 = scatter(thisx, fitmat(fitmat(:,3)==v & fitmat(:,6)==p,1),...
            'filled');
        h2.MarkerFaceColor = cm(v,:) + ( (1-cm(v,:))*0.5 );
        h2.SizeData = 8;        
        pp = pp + 1;        
        xlim([0 3]);
        title([voilist{v} ' ' parlabel{p}]);

        % plot best fit CV line
        hold all;
        thiscoef = bestcoefs_mean{p,v};
        if isempty(thiscoef)
            continue;
        end
        fplot = polyval(thiscoef,sort(thisx(1:10:end)));
        plot(sort(thisx(1:10:end)),fplot,'-','LineWidth',3,'Color',cm(v,:));        
        % plot the 95% CIs on the coefs
        for i = 1:2
            if size(bestcoefs_ci{p,v},1) == 1
                cdat = bestcoefs_ci{p,v}(i);
            else
                cdat = bestcoefs_ci{p,v}(i,:);
            end
            ftmp = polyval(cdat,sort(thisx(1:10:end)));
            plot(sort(thisx(1:10:end)),ftmp,'--','LineWidth',2,...
                'Color',cm(v,:));
        end

        clear thisx thiscoef fplot
    end
end

% make scatter plot of significant ones only

cm = lines(nv);
thisx = 0:0.25:2.5;

figure; pp = 1;
for p = 1:npar
    h(p) = subplot(1,npar,p);
    for v = 1:nv
        thiscoef = bestcoefs_mean{p,v};
        fplot = polyval(thiscoef,thisx);

        if isnan(bestmods(p,v))
            continue;
        elseif fdrmask(p,v)
            plot(thisx,fplot,'-','LineWidth',3,'Color',cm(v,:)); hold all;

            xlim([0 2.5]);
            title([parlabel{p}]);
        else
            plot(thisx,fplot,'LineStyle',':','LineWidth',2,'Color',cm(v,:)); hold all;

            xlim([0 2.5]);
            title([parlabel{p}]);
        end
    end
end