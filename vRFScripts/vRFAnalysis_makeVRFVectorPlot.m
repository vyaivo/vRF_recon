% vRFanalysis_makeVRFVectorPlot.m
% VAV 9/16/2015, cleaned for OSF 12/12/2016

% changelog - VAV 1/2/2017 edited to remove diff score outliers

clear xcomp ycomp fixs fixr;

root = load_root;
fitdir = 'vRFfits';
fn = fullfile(root,fitdir,'AllSubs_CompiledvRFs.mat');
load(fn);
condlist = {'left','right'};

% also want to load outliers to pull those out.
load(sprintf('%s%s/all_vRF_diffscores_noOutliers.mat',root,fitdir),'allout');

%%
% divide the vRFs into ecc bins
rstep = 0.3;
tharray = -pi:pi/4:pi;
ntbins = length(tharray)-1;
rarray = 0:rstep:5;
nrbins = length(rarray)-1;
nbins = ntbins*nrbins;

% make array of all bins
[t, r] = meshgrid(tharray(1:end-1),rarray(1:end-1));
% also do this but for bin centers
tcenter = tharray(1:end-1) + diff(tharray)/2;
rcenter = rarray(1:end-1) + diff(rarray)/2;

[tlab, rlab] = meshgrid(tcenter, rcenter);

% the imagesc/reshape process always generated upside down images because
% the y axis is flipped. so up until now everything has been
% kept upside down...including fitting to the image. so now here let's flip
% because the attended location was actually slightly below the midline
atlocs(:,2) = atlocs(:,2)*-1;

%%
for s = 1:length(sublist)
    for v = 1:length(voilist)
        newdat = alldat{s,v};
        newdat(allout{s,v},:,:) = [];   % remove outliers
        
        newdat(:,:,2) = newdat(:,:,2) * -1;
        
        xx = newdat(:,3,1);
        yy = newdat(:,3,2);
        thisdat = nan(size(newdat,1),2,size(newdat,3)-1);
        for c = 1:2
            theta = atan2( yy-atlocs(c,2), xx-atlocs(c,1) );
            rho = sqrt( (xx-atlocs(c,1)).^2 + (yy-atlocs(c,2)).^2 );
            % save the fixation distance from attend data
            fixdist(c,:) = rho;
            % save the distance from attend location data
            xc = newdat(:,c,1);
            yc = newdat(:,c,2);
            thisr = sqrt( (xc-atlocs(c,1)).^2 + ...
                (yc-atlocs(c,2)).^2 );
            % save dist param: attenddist - fixdist
            thisdat(:,c,1) = thisr - rho;
            % now sort into bins. also find the avg x/y vector length so we
            % can make our vec plots
            [~,~,tt] = histcounts(theta,tharray);
            [~,~,rr] = histcounts(rho,rarray);
            bininds = zeros(length(tt),2);
            for ti = 1:ntbins
                for ri = 1:nrbins
                    inds = find(tt == ti & rr == ri);
                    % bininds col 1: theta; col 2: rho
                    bininds(inds,1) = ti;
                    bininds(inds,2) = ri;
                    if isempty(inds)
                        xcomp(s,v,c,ti,ri) = NaN;
                        ycomp(s,v,c,ti,ri) = NaN;
                    else
                        xcomp(s,v,c,ti,ri) = mean(xc(inds) - xx(inds));
                        ycomp(s,v,c,ti,ri) = mean(yc(inds) - yy(inds));
                        alls(s,v,c,ti,ri) = mean(newdat(inds,c,3));
                        allr(s,v,c,ti,ri) = mean(newdat(inds,c,4)) - mean(newdat(inds,c,5));
                    end
                    clear inds
                end
            end
            
            for ti = 1:ntbins
                for ri = 1:nrbins
                    inds = find(tt == ti & rr == ri);                    
                    fixs(s,v,c,ti,ri) = mean(newdat(inds,3,3));
                    fixr(s,v,c,ti,ri) = mean(newdat(inds,3,4)) - mean(newdat(inds,3,5));
                end
            end
            
            clear xc yc bininds fixdist
        end
        
        clear newdat
    end
end

%% collapse across condition & plot

figure;

% for v = 1:length(voilist)
for v = 5;
    xcomps = squeeze(nanmean(nanmean(xcomp(:,v,:,:,:),1),3))';
    ycomps = squeeze(nanmean(nanmean(ycomp(:,v,:,:,:),1),3))';
    
%     ax = subplot(2,3,v);
    ax = gca;
    % set up data for vector plot
    allorigs = [(rlab(:).*cos(tlab(:)))+atlocs(1,1), ...
        (rlab(:).*sin(tlab(:)))+atlocs(1,2)]';
    h = quiverc_away_toward(allorigs(1,:),allorigs(2,:),...
        xcomps(:)',ycomps(:)',atlocs(1,:));
    hold all;
    plot(0,0,'k.');
    plot(atlocs(1,1),atlocs(1,2),'ko','LineWidth',2);
    plot(atlocs(2,1),atlocs(2,2),'o','LineWidth',2,'MarkerEdgeColor',[0.4 0.4 0.4]);                
    axis equal; box off;
    xlim([-5 5]); ylim([-4 4]);

    ax.FontSize = 10;
    title(voilist{v});
end

%% make plot of all ROIs / conditions
figure;
pp = 1;

for v = 1:length(voilist)
    for c = 1:2
        ax = subplot(3,4,pp);
        if c == 1
            ic = 2;
        else
            ic = 1;
        end

        xcomps = squeeze(nanmean(xcomp(:,v,c,:,:),1))';
        ycomps = squeeze(nanmean(ycomp(:,v,c,:,:),1))';
        % set up data for vector plot
        allorigs = [(rlab(:).*cos(tlab(:)))+atlocs(c,1), ...
            (rlab(:).*sin(tlab(:)))+atlocs(c,2)]';
        h = quiverc_away_toward(allorigs(1,:),allorigs(2,:),...
            xcomps(:)',ycomps(:)',atlocs(c,:));
        hold all;
        plot(0,0,'k.');
        plot(atlocs(c,1),atlocs(c,2),'ko','LineWidth',2);
        plot(atlocs(ic,1),atlocs(ic,2),'o','LineWidth',2,'MarkerEdgeColor',[0.4 0.4 0.4]);                
        axis equal; box off;
        xlim([-5 5]); ylim([-4 4]);
        
        ax.FontSize = 10;
        pp = pp+1;

        title([voilist{v} ', attend ' condlist{c}]);
    end
    clear cx cy
end