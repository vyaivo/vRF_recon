%%% vRFanalysis_attnModulations_stats
%%% VAV 6/9/2016, edited for OSF 12/12/2016

%%% fits a flat line, sloped line, and a polynomial to the vRF x
%%% attend_dist relationships. Then performs a nested F-test to see which
%%% is the best descriptor. Bootstraps the best fit parameters to get 95%
%%% CIs on all of them.

clear all;

%% load the data
root = load_root;
fitdir = 'vRFfits';
fn = fullfile(root,fitdir,'AllSubs_CompiledvRFs.mat');
load(fn);
savefn = fullfile(root,fitdir,'vRF_attnMods_hemi_stats.mat');

parlabel = {'position','size','amplitude','baseline'};

%% sort the vRFs

% to run the ANOVA, you need giant vectors with labels. so just slightly
% shift the data around to achieve that...

% dist from attend bins
distbins = 0:0.25:2.5;

adat = [];slab = [];vlab = [];dlab = [];blab = [];
for s = 1:ns
    for v = 1:nv
        if isempty(alldat{s,v})
            continue;
        end
        
        thisdat = alldat{s,v}(:,:,1:npar);
        
        %% calculate vRF param changes (attend - fix)
        
        attndat = nan(size(thisdat,1),2,size(thisdat,3)-1);
        for c = 1:nc-1
            
            % distance from attention locus during fix cond
            dist_at_fix(c,:) = sqrt( (thisdat(:,3,1)-atlocs(c,1)).^2 + ...
                (thisdat(:,3,2)-atlocs(c,2)).^2 );            
            % sort vRFs by distance from attention locus during fix cond
            [~,bi(c,:)] = histc(dist_at_fix(c,:),distbins);
            
            % distance from attention locus during attention cond
            dist_attend = sqrt( (thisdat(:,c,1)-atlocs(c,1)).^2 + ...
                (thisdat(:,c,2)-atlocs(c,2)).^2 );
            
            % now calc diff between position measures above
            attndat(:,c,1) = dist_attend' - dist_at_fix(c,:);
            
            for p = 3:npar
                % calculate difference from fixation condition
                attndat(:,c,p-1) = thisdat(:,c,p) - thisdat(:,3,p);
            end
            
            clear dist_attend
        end
        
        % we want to fold the data together so we should get 2x the vRFs
        catdat = squeeze(cat(1,attndat(:,1,:),attndat(:,2,:)));
        catbin = [bi(1,:) bi(2,:)];
        catfdist = [dist_at_fix(1,:), dist_at_fix(2,:)];
        
        % thisdat has size nvRF x npar
        sl = repmat(s,size(catdat));                % subject label
        vl = repmat(v,size(catdat));                % VOI label
        bl = repmat(catbin',1,size(catdat,2));
        dl = repmat(catfdist',1,size(catdat,2));    % distance from attend (continuous covariate)
        
        % save out all the labels for the ANOVA
        slab = cat(1,slab,sl);
        vlab = cat(1,vlab,vl);
        blab = cat(1,blab,bl);
        dlab = cat(1,dlab,dl);
        adat = cat(1,adat,catdat);
                
        clearvars dist_at_fix catdat catbin catfdist  
        clearvars newdat attndat bi sl vl bl dl
    end
end

%% now make the big fit matrix

% update number of conditions & params now that we've fiddled a bit
nc = 1; npar = npar-1;

plab = nan(size(adat));
for p = 1:npar
    plab(:,p) = p;
end
fitmat = [adat(:) slab(:) vlab(:) blab(:) dlab(:) plab(:)];
% fitmat is: data, sub label, VOI label, bin label, distance covariate, param label

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% GET MEAN VRF CHANGES
rs2 = sum(clock*10000);
rng(rs2);
niters = 10000;

all_diffscores = nan(npar,nv,niters);
mean_diffscores = nan(npar,nv);
ci_diffscores = nan(npar,nv,2);
pval_diffscores = nan(npar,nv);

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

%%
for p = 1:npar
    for v = 1:nv
        
        %% first fit each model using 50/50 cross-validation
        thisd = fitmat(fitmat(:,3)==v & fitmat(:,6)==p,:);
        
        db = thisd(:,5);
        dat = thisd(:,1);
        nanb = isnan(dat);
        dat(nanb) = [];
        db(nanb) = [];
        
        ssres = nan(niters,3); sstot = nan(niters,3); r2 = nan(niters,3);
        trncv = logical(zeros(numel(dat),niters)); tstcv = trncv;
        parfor i = 1:niters
            [trncv(:,i),tstcv(:,i)] = crossvalind('HoldOut',numel(dat),0.5);
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

%%
coef_pval = nan(npar,nv);
bestcoefs_mean = cell(npar,nv);
bestcoefs_ci = cell(npar,nv);

for p = 1:npar
    for v = 1:nv
        
        bestmod = max(find(squeeze(allpvals(p,v,:)) < .05));
        if isempty(bestmod)
            continue;
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

%% save

save(savefn,'rs','allssres','allr2s','bestmods','bestr2s','bestcoefs_mean',...
    'bestcoefs_ci','allFvals','allpvals','niters','coefs','coef_pval',...
    'fdrthresh','fdrmask','rs2','all_diffscores','mean_diffscores',...
    'ci_diffscores','pval_diffscores','-v7.3');

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
        fplot = polyval(thiscoef,thisx);
        plot(thisx,fplot,'-','LineWidth',3,'Color',cm(v,:));        
        
        clear thisx thiscoef fplot
    end
end

%% make scatter plot of significant ones only

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
