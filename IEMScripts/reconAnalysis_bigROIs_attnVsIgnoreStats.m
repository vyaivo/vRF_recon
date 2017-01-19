%%% reconAnalysis_bigROIs_attnVsIgnoreStats.m
%%% does randomization stats on fits to reconstructions of all 51 positions
%%% of the stimulus. this version bins the positions by eccentricity from
%%% the attend location.

%%% VAV 4/16/2015
%%% cleaned for OSF 12/14/2016

% dir stuff
root = load_root;
rdir = 'recons';
fitdir = 'reconfits';
fitfix = 'acrossSess'; trnfix = 'TrnAllJitter';
sublist = {'AA','AI','AL','AP','AR','AT','AU'};
voilist = {'superVis','superParietal'};
% voilist = {'V1','V2','V3','V3AB','V4','IPS0'};
ns = length(sublist);
nv = length(voilist);

% set up some n
niters = 10000;

statfn = fullfile(root,fitdir,'reconAnalysis_bigROIs_statswANOVA.mat');

%%

datfn = sprintf('%s%s/AllSubs_CompiledRecons_bigROIs_%s.mat',root,fitdir,trnfix);
load(datfn);
alls = size(reconFits);
clearvars allRecons reconSubAvg      % these are images; we don't need them

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% calculate mean recon position err %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% e.g., how well do the stimulus reconstructions predict the stim pos?
% accompanies Fig 4B citation in text.

for pp = 1:size(locs,1)
    x = squeeze(reconFits(:,:,:,pp,1));
    y = squeeze(reconFits(:,:,:,pp,2));
    reconerr(:,:,:,pp) = sqrt( (x-locs(pp,1)).^2 + (y-locs(pp,1)).^2 );
end

rs1 = 496740;
rng(rs1);
% now resample so we have CIs on pos errs
reserr = nan(niters,nv,nc,size(locs,1));
parfor i = 1:niters
    rsi = randi(ns,ns,1);
    reserr(i,:,:,:) = mean(reconerr(:,rsi,:,:),2);
end

mean_recerr = squeeze(mean(reserr,1));
ci_recerr = prctile(reserr,[2.5 97.5]);

mean_overallerr = mean(reshape(reserr(:,1:nv,:,:),1,[]));
ci_overallerr = prctile(reshape(reserr(:,1:nv,:,:),1,[]),[2.5 97.5]);

save(statfn,'reconerr','rs1','mean_overallerr','ci_overallerr','-v7.3');

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% run 2 way ANOVA of dist fr attn loc x recon fit param %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% change x & y --> center dist from attended pos

alls2 = size(rdat);
npar = alls2(5)-1;
rdat2 = nan([alls2(1:4),alls2(5)-1]);
rdat2(:,:,:,:,2:4) = rdat(:,:,:,:,3:5);
for b = 1:nbins
    lbi = find(eccbini(1,:) == b);
    rbi = find(eccbini(2,:) == b);
    % first collapse the x & y centers into dist from attend/ignore locs
    atleft = sqrt( (reconFits(:,:,1,lbi,1)-attendlocs(1,1)).^2 + ...
        (reconFits(:,:,1,lbi,2)-attendlocs(1,2)).^2 );   % dist from attend left
    atright = sqrt( (reconFits(:,:,2,rbi,1)-attendlocs(2,1)).^2 + ...
        (reconFits(:,:,2,rbi,2)-attendlocs(2,2)).^2 );   % dist from attend right
    atcenters = mean(cat(4,atleft,atright),4);
    
    % the data from attend left runs, sorted from the distance from the attend
    % right locus
    igleft = sqrt( (reconFits(:,:,1,rbi,1)-attendlocs(2,1)).^2 + ...
        (reconFits(:,:,1,rbi,2)-attendlocs(2,2)).^2 );
    % the data from attend right runs, sorted from the distance from the attend
    % left locus
    igright = sqrt( (reconFits(:,:,2,lbi,1)-attendlocs(1,1)).^2 + ...
        (reconFits(:,:,2,lbi,2)-attendlocs(1,2)).^2 );   % dist from attend left
    igcenters = mean(cat(4,igleft,igright),4);
    % now save these to rdat2
    rdat2(:,:,1,b,1) = atcenters;
    rdat2(:,:,2,b,1) = igcenters;
    % take the mean of these positions near the attend locations
    rdat2(:,:,1,b,2:npar) = mean(cat(4,reconFits(:,:,1,lbi,3:npar+1),reconFits(:,:,2,rbi,3:npar+1)),4);
    % take the mean of all these positions near the ignored location
    rdat2(:,:,2,b,2:npar) = mean(cat(4,reconFits(:,:,1,rbi,3:npar+1),reconFits(:,:,2,lbi,3:npar+1)),4);
    clear lbi rbi
end

%% now make a big data matrix w/ labels
% this is needed for the ANOVA

% col 1: fit data
fitvec = rdat2(:);
% col 2: VOI labels
vlabs = nan(size(rdat2));
for v = 1:nv
    vlabs(v,:,:,:,:) = v;
end
vvec = vlabs(:);
% col 3: sub labels
slabs = nan(size(rdat2));
for s = 1:ns
    slabs(:,s,:,:,:) = s;
end
svec = slabs(:);
% col 4: attend labels
alabs = nan(size(rdat2));
for a = 1:2
    alabs(:,:,a,:,:) = a;
end
avec = alabs(:);
% col 5: bin labels
blabs = nan(size(rdat2));
for b = 1:nbins
    blabs(:,:,:,b,:) = b;
end
bvec = blabs(:);
% col 6: fit param
plabs = nan(size(rdat2));
for p = 1:npar
    plabs(:,:,:,:,p) = p;
end
pvec = plabs(:);
% now concatenate everything!
fitmat = [fitvec vvec svec avec bvec pvec];

clear vlabs slabs clabs blabs plabs

%% OMNIBUS TEST

rs2 = 6078943;
rng(rs2);
pom = nan(npar,nv);
f_om_shuf = nan(npar,nv,niters);
f_om_real = nan(npar,nv,niters);

for p = 1:npar
    for v = 1:nv
        % only get the data where fit param = p & VOI = v
        pvinds = (fitmat(:,2) == v) & (fitmat(:,6) == p);
        thisd = fitmat(pvinds,:);
        % collapse condition & posbin
        cp = [thisd(:,4) thisd(:,5)];
        [~,~,newlabs] = unique(cp,'rows');
        realf = RMAOV1([thisd(:,1), newlabs thisd(:,3)],0.05);
        % DO THE SHUFFLE
        tic
        parfor i = 1:niters
            randmat = nan(size(thisd,1),2);
            nts = size(randmat,1)/ns;
            % randomize WITHIN subject
            for s = 1:ns
                is = find(thisd(:,3) == s);
                sp = randperm(length(is));
                randmat((s-1)*nts+1:(s-1)*nts+nts,:) = [newlabs(is(sp)),...
                    thisd(is(sp),3)];
            end
            [fom(i), ~] = RMAOV1([thisd(:,1), randmat],0.05);
        end
        toc
        pom(p,v) = sum(realf <= fom) / niters;
        if pom(p,v) == 0
            pom(p,v) = 1/(niters*10);
        end
        f_om_shuf(p,v,:) = fom;
        f_om_real(p,v,:) = realf;
        
        clear thisd cp newlabs fom realf
    end
    % now calc the fdr correction
    [pfdr(p,:), pmasked(p,:)] = fdr(squeeze(pom(p,:)), 0.05);
end

%%

save(statfn,'reconerr','rs1','mean_overallerr','ci_overallerr','fitmat',...
    'f_om_real','f_om_shuf','binlabels','eccbini','pfdr','pom','pmasked',...
    'rs2','npar','nv','niters','-v7.3');

%% ALL F TESTS

rs3 = 968056;
rng(rs3);

for p = 1:npar
    for v = 1:nv
        % only get the data where fit param = p & VOI = v
        pvinds = (fitmat(:,2) == v) & (fitmat(:,6) == p);
        thisd = fitmat(pvinds,:);
        % collapse condition & posbin
        ftab = rm_anova2(thisd(:,1),thisd(:,3),thisd(:,4),thisd(:,5),...
            {'cond','posbin'});
        realf_cond = ftab{2,5};
        realf_pos = ftab{3,5};
        realf_condbypos = ftab{4,5};
        % DO THE SHUFFLE
        f_cond_shuf = nan(niters,1);
        f_pos_shuf = nan(niters,1);
        f_condbypos_shuf = nan(niters,1);
        tic
        parfor i = 1:niters
            rf1 = nan(size(thisd,1),1);
            rf2 = rf1;
            rs = rf1;
            nts = size(rf1,1)/ns;
            % randomize WITHIN subject
            for s = 1:ns
                is = find(thisd(:,3) == s);
                sp = randperm(length(is));
                rf1((s-1)*nts+1:(s-1)*nts+nts,:) = thisd(is(sp),4);
                rf2((s-1)*nts+1:(s-1)*nts+nts,:) = thisd(is(sp),5);
                rs((s-1)*nts+1:(s-1)*nts+nts,:) = thisd(is(sp),3);
            end
            ftabr = rm_anova2(thisd(:,1),rs,rf1,rf2,{'cond','posbin'});
            f_cond_shuf(i) = ftabr{2,5};
            f_pos_shuf(i) = ftabr{3,5};
            f_condbypos_shuf(i) = ftabr{4,5};
        end
        toc
        pcond(p,v) = sum(realf_cond <= f_cond_shuf) / niters;
        ppos(p,v) = sum(realf_pos <= f_pos_shuf) / niters;
        pcondbypos(p,v) = sum(realf_condbypos <= f_condbypos_shuf) / niters;
        if pcond(p,v) == 0
            pcond(p,v) = 1/(niters*10);
        end
        if ppos(p,v) == 0
            ppos(p,v) = 1/(niters*10);
        end
        if pcondbypos(p,v) == 0
            pcondbypos(p,v) = 1/(niters*10);
        end
        all_f_cond_shuf(p,v,:) = f_cond_shuf;
        all_f_condbypos_shuf(p,v,:) = f_condbypos_shuf;
        all_f_pos_shuf(p,v,:) = f_pos_shuf;
        all_f_cond(p,v,:) = realf_cond;
        all_f_pos(p,v,:) = realf_pos;
        all_f_condbypos(p,v,:) = realf_condbypos;
        
        clear pvinds thisd cp newlabs realf_cond realf_pos realf_condbypos
        clear f_cond_shuf f_condbypos_shuf
    end
end

% calc fdr correction
for p = 1:npar
    [pfdr_cond(p), pmask_cond(p,:)] = fdr(pcond(p,:), 0.05);
    [pfdr_pos(p), pmask_pos(p,:)] = fdr(ppos(p,:), 0.05);
    [pfdr_condbypos(p), pmask_condbypos(p,:)] = fdr(pcondbypos(p,:), 0.05);
end

%%
save(statfn,'reconerr','rs1','mean_overallerr','ci_overallerr','fitmat',...
    'f_om_real','f_om_shuf','binlabels','eccbini','pfdr','pom','pmasked',...
    'rs2','npar','nv','niters','all_f_cond','all_f_pos','all_f_condbypos','rs3',...
    'all_f_cond_shuf','all_f_pos_shuf','all_f_condbypos_shuf',...
    'pfdr_cond','pfdr_pos','pfdr_condbypos','pmask_cond',...
    'pmask_pos','pmask_condbypos','pcond','ppos','pcondbypos','-v7.3');

%% bootstrap the bin data so we get real means / CIs
rs4 = 213096;
rng(rs4);

all_iterbins = nan(npar,nv,nc,niters,nbins);
mean_reconbin = nan(npar,nv,nc,nbins);
ci_reconbin = nan(npar,nv,nc,nbins,2);

for p = 1:npar
    for v = 1:nv
        for c = 1:nc
            % only get the data where fit param = p & VOI = v
            pvinds = (fitmat(:,2) == v) & (fitmat(:,6) == p) & ...
                (fitmat(:,4) == c);
            thisd = fitmat(pvinds,:);
            iterbin = nan(niters,nbins);

            tic
            parfor i = 1:niters            
                rd = nan(size(thisd,1),1);
                rb = rd;
                nts = size(rd,1)/ns;
                % resample WITHIN subject
                for s = 1:ns
                    is = find(thisd(:,3) == s);
                    sp = randi(length(is),length(is),1);
                    rd((s-1)*nts+1:(s-1)*nts+nts,:) = thisd(is(sp),1);
                    rb((s-1)*nts+1:(s-1)*nts+nts,:) = thisd(is(sp),5);
                end
                for b = 1:nbins
                    iterbin(i,b) = nanmean(rd(rb==b));
                end
            end
            toc

            all_iterbins(p,v,c,:,:) = iterbin;
            mean_reconbin(p,v,c,:) = nanmean(iterbin,1);
            ci_reconbin(p,v,c,:,:) = prctile(iterbin,[2.5,97.5])';
        end
    end
end

%%
save(statfn,'reconerr','rs1','mean_overallerr','ci_overallerr','fitmat',...
    'f_om_real','f_om_shuf','binlabels','eccbini','pfdr','pom','pmasked',...
    'rs2','npar','nv','niters','all_f_cond','all_f_pos','all_f_condbypos',...
    'rs3','all_f_cond_shuf','all_f_pos_shuf','all_f_condbypos_shuf',...
    'pfdr_cond','pfdr_pos','pfdr_condbypos','pmask_cond',...
    'pmask_pos','pmask_condbypos','pcond','ppos','pcondbypos','rs4',...
    'all_iterbins','mean_reconbin','ci_reconbin','-v7.3');

%% get some stuff ready to plot
% Fig 5

if ~exist('fitmat','var')
    load(statfn);
end
nbins = length(binlabels);
nc = 2;

yparlab = {'distance from attend location','size','amplitude','baseline'};
cp = lines;
ccol = [cp(1,:); cp(2,:)];

%% actually plot
pp = 1;
figure;
for p = 1:npar
    for v = 1:nv
        ah(v) = subplot(npar,nv,pp); hold all;
        for c = 1:nc
            if p == 2
                datc = rad2fwhm(squeeze(mean_reconbin(p,v,c,:)));
                errdatc = abs(rad2fwhm(squeeze(ci_reconbin(p,v,c,:,:))) ...
                    - datc);
            else
                datc = squeeze(mean_reconbin(p,v,c,:));
                errdatc = abs(squeeze(ci_reconbin(p,v,c,:,:)) - datc);
            end
                
            h = errorbar(1:nbins,datc,errdatc(:,1),errdatc(:,2),...
                'Color',ccol(c,:));
            h.LineWidth = 2;
            ax = gca;
            ax.XTick = 1:nbins;
            ax.XTickLabel = binlabels;        
            ax.XTickLabelRotation = -65;
            ax.LineWidth = 1;
            ax.FontSize = 10;
            xlim([0 nbins+1]);
        end
        match_ylim(ah);
        yd = ylim;
        if pmask_cond(p,v); plot(2, yd(2)*0.9, 'k*'); end;
        if pmask_pos(p,v); plot(3, yd(2)*0.9, 'ko'); end;
        if pmask_condbypos(p,v); plot(4, yd(2)*0.9, 'kx'); end;
        if v == 1
            ylabel(yparlab{p});
        end
        
        if p == 1
            title(voilist{v});
        end
        pp = pp + 1;
    end
    match_ylim(ah); clear ah yd;
end

legend('attend','ignore');