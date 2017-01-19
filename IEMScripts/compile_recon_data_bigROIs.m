function compile_recon_data_bigROIs(trnfix)
% compile_recon_data_allsubs.m
% VAV 4/16/2015
% puts all recon data together in a big matrix.
% cleaned up for OSF 12/14/2016

if nargin < 1
    trnfix = 'TrnAllJitter';
end

%% dir stuff
root = load_root;
rdir = 'recons';
fitdir = 'reconfits';
fitfix = 'acrossSess';
sublist = {'AA','AI','AL','AP','AR','AT','AU'};
ns = length(sublist);
voilist = {'superVis','superParietal'};
nv = length(voilist);
% voilist = {'V1','V2','V3','V3AB','V4','IPS0'};
% nv = length(voilist);

savefn = sprintf('%s%s/AllSubs_CompiledRecons_bigROIs_%s.mat',root,fitdir,trnfix);
attendlocs = [-2.1074, 0.3664; 2.1233, 0.3664];
npos = 51;
npar = 5;
nc = 2;
skipsub = cell(nv,1);

%% now load each file & stick it in the big matrix

reconFits = nan(nv,ns,nc,npos,npar);

for v = 1:length(voilist)
    for s = 1:length(sublist)
        %% load recon file
        fn1 = sprintf('%s%s/%s_%s_%s_%s.mat',root,rdir,sublist{s},...
            trnfix,fitfix,voilist{v});
        if ~exist(fn1,'file')
            fprintf('Missing file: %s\n',fn1);
            skipsub{v} = cat(1,skipsub{v},s);
            continue;
        end
        load(fn1,'recon_avg','locs','res','fov');
        allRecons(s,v,:,:,:) = recon_avg;
        
        %% load fit file
        fn2 = sprintf('%s/%s/%s_%s_ReconGridFits_%s_%s.mat',root,fitdir,...
            sublist{s},voilist{v},fitfix,trnfix);
        if ~exist(fn2,'file')
            fprintf('Sorry, cannot find %s...\n',fn2);
            continue;
        end
        load(fn2);
        
        reconFits(v,s,:,:,:) = bfpar;
        fitErr(v,s,:,:) = bferr;        
        
        clear bferr bfpar
    end
end

reconSubAvg = squeeze(nanmean(allRecons));

%% sort into hemifield & eccentricity bins

% 1) define the bins
% We will bin the data by distance from the attended location. Since the 
% stimulus positions form a regular triangular grid, with one point being
% the attended location, we will get roughly iso-eccentric rings around the
% stimulus. These bin definitions are created to respect the iso-eccentric
% boundaries and then to create roughly equal numbers of reconstructions in
% each bin (once those boundaries break down).
binedges = [0 0.1 0.86 1.48 1.70 2.25 2.55 3.1];    % 0 & 3.1 are catch bins to get rid of any recons that fall outside our other bins
binlabels = [0 0.85 1.47 1.69 2.24 2.54];           % these are the actual centers of the bins
nbins = length(binedges)-2;

% 2) calculate distance from attended location, and then sort into bins
rdat = nan([nv,ns,nc,nbins,npar]);
locdist = nan(2,npos);
for c = 1:2
    % calculate dist from attend location
    locdist(c,:) = sqrt( (locs(:,1)-attendlocs(c,1)).^2 + ...
        (locs(:,2)-attendlocs(c,2)).^2 );
    % bin all positions
    [nb(c,:),eccbini(c,:)] = histc(locdist(c,:),binedges);
end

% 3) now collapse over attend L/R
% now loop over bins and combine into attended / unattended hemifield

% don't take the bins at the end (they're trash -- just other locations
% that fall outside 2.55, which starts crossing into the other
% hemifield)

for b = 1:nbins
    lbi = find(eccbini(1,:) == b);
    rbi = find(eccbini(2,:) == b);
    % take the mean of these positions near the attend locations
    rdat(:,:,1,b,:) = mean(cat(4,reconFits(:,:,1,lbi,:),reconFits(:,:,2,rbi,:)),4);
    % take the mean of all these positions near the ignored location
    rdat(:,:,2,b,:) = mean(cat(4,reconFits(:,:,1,rbi,:),reconFits(:,:,2,lbi,:)),4);
end

%% resample fit errors
rs = 56904;
rng(rs);
niters = 10000;
err_dist = nan(1,niters);
parfor i = 1:niters
    rsi = randsample(numel(fitErr),numel(fitErr),1);
    err_dist(i) = nanmean(fitErr(rsi));
end
mean_fitErr = mean(err_dist);
ci_fitErr = prctile(err_dist,[2.5,97.5]);

%%

save(savefn,'allRecons','reconSubAvg','reconFits','rdat','binlabels','eccbini',...
    'sublist','voilist','attendlocs','locdist','locs','nbins','ns','nv',...
    'nbins','nc','npar','fov','res','fitErr','skipsub','-v7.3');
fprintf('Saved %s!\n',savefn);