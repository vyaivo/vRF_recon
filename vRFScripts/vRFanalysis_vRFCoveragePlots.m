%%% vRFanalysis_vRFCoveragePlots.m
%%% VAV 12/15/2015

root = load_root;
fitdir = 'vRFfits';

%% load the data
fn = fullfile(root,fitdir,'AllSubs_CompiledvRFs.mat');
load(fn,'alldat','allerr','allgoodvox','sublist','voilist','atlocs',...
    'npar','nc','ns','nv','vrf_threshhold_table');
load(sprintf('%s%s/all_vRF_diffscores_noOutliers.mat',root,fitdir),'allout');

%%
x = linspace(-6,6,101);
y = linspace(-6,6,101);
[xx,yy] = meshgrid(x,y);

%%
climits = [0 40];
% now make some plots...
pp=1; clear h;
nodat = [];
for s = 1:ns
    for v = 1:nv        
        if isempty(alldat{s,v})
            nodat = [nodat pp];
            pp = pp+1;
            continue;
        end
        % remove outliers
        alldat{s,v}(allout{s,v},:,:) = [];
        
        h(pp) = subplot(ns,nv,pp); hold all;
        mux = squeeze(alldat{s,v}(:,3,1));
        muy = squeeze(alldat{s,v}(:,3,2))*-1;
        % convert size FWHM to size constant
        sz = alldat{s,v}(:,3,3) ./ rad2fwhm(1);
        for vv = 1:length(mux)
            dat(vv,:,:) = flipud(make2dcos(xx,yy,mux(vv),muy(vv),...
                sz(vv),7));
        end
        sdat = squeeze(nansum(dat,1));
        prfcoverage(s,v,:,:) = sdat;
        
        colormap parula;
        imagesc(x,-1*y,sdat);
        sh = scatter(mux, muy,'.');
        sh.CData = [0.5 0.5 0.5];
        sh.SizeData = 12;
        axis equal; axis off;
        xlim([-6 6]); ylim([-6 6]);
        
        if pp <= nv
            title(voilist{v});
        end
        
        clear mux muy sdat dat
        pp = pp+1;
    end
end
h(nodat) = [];
ch = match_clim(h,climits);

%% all subjects
allx = []; ally = [];
for s = 1:ns
    for v = 1
        allx = cat(1,allx,squeeze(alldat{s,v}(:,3,1)));
        ally = cat(1,ally,squeeze(alldat{s,v}(:,3,2)));
    end
end

allcov = squeeze(sum(prfcoverage(:,1,:,:)));
figure;
imagesc(x,y,allcov,ch); colorbar;
hold on; axis equal; box off; axis off;
sh2 = scatter(allx,ally);
sh2.CData = [0.5 0.5 0.5];