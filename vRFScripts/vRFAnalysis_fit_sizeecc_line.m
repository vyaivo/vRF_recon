% vRFanalysis_fit_sizeecc_line.m
% VAV cleaned for OSF 12/12/2016
% changelog: VAV 1/2/2017 removed diff score outliers

%% path stuff
root = load_root;
fitdir = 'vRFfits';
fn = fullfile(root,fitdir,'AllSubs_CompiledvRFs.mat');
load(fn);
load(sprintf('%s%s/all_vRF_diffscores_noOutliers.mat',root,fitdir),'allout');
savefn = sprintf('%s%s/%s',root,fitdir,'vrf_size_ecc_data.mat');

%% analysis params

niters = 10000;              % bootstrap iterations
rs = sum(clock*10000);
rng(rs);
eccbins = 0:0.25:2.25;
mbins = eccbins(1:end-1) + (diff(eccbins)/2);
X = [mbins; ones(size(mbins))]';
ycoef_boot = nan(nv,niters,2);
yline_boot = nan(nv,niters,numel(mbins));
ysz_boot = nan(nv,niters,numel(mbins));
yecc_boot = nan(nv,niters,numel(mbins));

%% do the regression with bootstrapping so we can get error bars
for v = 1:nv
    szdat = []; eccdat = []; slab = [];
    for s = 1:ns
        alldat{s,v}(allout{s,v},:,:) = [];
        szdat = cat(1,szdat,alldat{s,v}(:,3,3));
        eccdat = cat(1,eccdat,sqrt(alldat{s,v}(:,3,1).^2 + alldat{s,v}(:,3,2).^2));
        slab = cat(1,slab,repmat(s,size(alldat{s,v},1),1));
    end
    nt = numel(szdat);
    for i = 1:niters
        ri = randsample(nt,nt,1);
        % resample within subjects
        sd = nan(size(szdat));
        ed = nan(size(eccdat));
        for s = 1:ns
            sl = find(slab == s);
            ri = randsample(numel(sl),numel(sl),1);
            sd(sl) = szdat(sl(ri));
            ed(sl) = eccdat(sl(ri));
        end
        [~,~,bi] = histcounts(ed,eccbins);
        sz = nan(size(mbins)); ec = sz;
        for ii = 1:max(bi)
            sz(ii) = nanmean(sd(bi==ii));
            ec(ii) = nanmean(ed(bi==ii));
        end
        ysz_boot(v,i,:) = sz;
        yecc_boot(v,i,:) = ec;
        ycoef_boot(v,i,:) = regress(sz',X);
        yline_boot(v,i,:) = mbins.*ycoef_boot(v,i,1) + ycoef_boot(v,i,2);
    end
    % now calc the 95% CI for sz, ecc, & fit
    sz_ci(v,:,:) = prctile(ysz_boot(v,:,:),[2.5 97.5]);
    ecc_ci(v,:,:) = prctile(yecc_boot(v,:,:),[2.5 97.5]);
    yline_ci(v,:,:) = prctile(yline_boot(v,:,:),[2.5 97.5]);
end

%% now also bootstrap the remaining bins...
eccbins2 = 2.25:0.25:5.5;
mbins2 = eccbins2(1:end-1) + (diff(eccbins2)/2);
y2sz_boot = nan(nv,niters,numel(mbins2));
y2ecc_boot = nan(nv,niters,numel(mbins2));
sz2_ci = nan(nv,2,numel(mbins2));
ecc2_ci = nan(nv,2,numel(mbins2));

for v = 1:nv
    szdat = []; eccdat = []; slab = [];
    for s = 1:ns
        szdat = cat(1,szdat,alldat{s,v}(:,3,3));
        eccdat = cat(1,eccdat,sqrt(alldat{s,v}(:,3,1).^2 + alldat{s,v}(:,3,2).^2));
        slab = cat(1,slab,repmat(s,size(alldat{s,v},1),1));
    end
    nt = numel(szdat);
    for i = 1:niters
        ri = randsample(nt,nt,1);
        
%         sd = szdat(ri);
%         ed = eccdat(ri);
%         [~,~,bi] = histcounts(ed,eccbins2);

        % resample within subjects
        sd = nan(size(szdat));
        ed = nan(size(eccdat));
        for s = 1:ns
            sl = find(slab == s);
            ri = randsample(numel(sl),numel(sl),1);
            sd(sl) = szdat(sl(ri));
            ed(sl) = eccdat(sl(ri));
        end
        [~,~,bi] = histcounts(ed,eccbins2);
        sz = nan(size(mbins2)); ec = sz;
        for ii = 1:max(bi)
            sz(ii) = nanmean(sd(bi==ii));
            ec(ii) = nanmean(ed(bi==ii));
        end
        y2sz_boot(v,i,:) = sz;
        y2ecc_boot(v,i,:) = ec;
    end
    % now calc the 95% CI for sz, ecc, & fit
    sz2_ci(v,:,:) = prctile(y2sz_boot(v,:,:),[2.5 97.5]);
    ecc2_ci(v,:,:) = prctile(y2ecc_boot(v,:,:),[2.5 97.5]);
    clear vind szdat eccdat nt
end


%% now test the significance of the fit slope
% first test if slope > 0
for v = 1:nv
    pslope(v) = 1 - mean(ycoef_boot(v,:,1)>0);
end

%% save

fprintf('SAVING...');
save(savefn,'ysz_boot','sz_ci','yline_boot','yline_ci',...
    'ysz_boot','sz_ci','y2sz_boot','sz2_ci','sublist','voilist','rs',...
    'mbins','mbins2','ycoef_boot','pslope','-v7.3');
fprintf('DONE!\n');

%% plot

% now plot that stuff.
cm = parula(nv+1);
figure;
for v = 1:nv
    szd = nanmean(squeeze(ysz_boot(v,:,:)));
    lowci = abs(squeeze(sz_ci(v,1,:))-szd');
    highci = abs(squeeze(sz_ci(v,2,:))-szd');
    h = errorbar(mbins,szd,lowci,highci,'o','MarkerFaceColor',cm(v,:),...
        'MarkerEdgeColor',cm(v,:),'Color',cm(v,:),'MarkerSize',8,...
        'LineWidth',2);
    hold on;
end
legend({'V1','V2','V3','V3A/B','V4','IPS0'});

for v = 1:nv
    yl = nanmean(squeeze(yline_boot(v,:,:)));
    lowy = abs(squeeze(yline_ci(v,1,:))-yl');
    highy = abs(squeeze(yline_ci(v,2,:))-yl');
    h = shadedErrorBar(mbins,yl,[highy,lowy]',{'LineWidth',2,'Color',cm(v,:)},1);
    h.LineWidth = 1;
    hold on;
end

box off;
xlabel('eccentricity (distance from fixation)');
ylabel('size (degrees)');

for v = 1:nv
    szd = nanmean(squeeze(y2sz_boot(v,:,:)));
    lowci = abs(squeeze(sz2_ci(v,1,:))-szd');
    highci = abs(squeeze(sz2_ci(v,2,:))-szd');
    h = errorbar(mbins2,szd',lowci,highci,'o--','MarkerFaceColor',cm(v,:),...
        'MarkerEdgeColor',cm(v,:),'Color',cm(v,:),'MarkerSize',8,...
        'LineWidth',2);
    hold on;
end
ax = gca; ax.FontSize = 9;