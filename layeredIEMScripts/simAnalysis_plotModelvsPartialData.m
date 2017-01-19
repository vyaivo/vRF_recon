% simAnalysis_plotModelvsPartialData.m
% VAV 6/16/2015
% cleaned for OSF 12/19/2016

root = load_root;
simdir = 'sims/';

% voilist = {'V1','V2','V3','V3AB','V4','IPS0'};
% nv = length(voilist);

nc = 2; np = 4; npos = 51;

% load layered IEM recon fits
bigfitfn = sprintf('%s%s/all_sim_covNoise_fitsWModel.mat',root,simdir);
load(bigfitfn);

for c = 1:nc
    % save out the actual dist from attend
    posdist(c,:) = sqrt( (locs(condlocs(c,:),1)-atlocs(c,1)).^2 + ...
        (locs(condlocs(c,:),2)-atlocs(c,2)).^2 );
end
% now bin the positions of the difference scores
binedges = [0 0.1 0.9 1.55 1.75 2.4 2.55 3.5];
binlabels = [0 0.85 1.47 1.69 2.24 2.54];
nbins = length(binedges)-2;
[nb,eccbini] = histc(posdist(1,:),binedges);

%%

% want to plot the real data as bar plots
% allmeanbdats is nvoi x ncond x nrecon x npar x nmodels

for v = 1:nv
    for p = 1:np
        for b = 1:nbins
            
            % get bin means for real data
            fulldat = reducedrec(v,:,eccbini==b,p);
            meanfullpar(v,p,b) = nanmean(fulldat(:));
            
            clear fulldat
        end
    end
end

%%
parlabels = {'position','size','amplitude','baseline'};
plotm = [1 2 3 4];
% plotm = [8 9 10];
% colors
cml = lines;
cml = cat(1,[0 0 0], cml);
% y limits
yl = [-1 1; -1 1; -0.5 0.75; -0.2 0.1];

figure;
pp = 1;
for p = 1:np
    for v = 1:nv
        h(v) = subplot(np,nv,pp); hold all;
        
        % plot real data
        hb = bar(1:nbins,squeeze(meanfullpar(v,p,:)));
        hb.FaceColor = [0.5 0.5 0.5];        
        
        % plot no sz + no shift models
        mi = 1;
        for m = plotm
            % mean over the bins
            mfun = @(xi) nanmean( reshape(modeldat(m,v,:,:,xi,p),1,[]) );
            mdat = table2array(rowfun(mfun,table(bsxfun(@eq,(1:nbins)',eccbini))));
%             hp = errorbar(1:nbins,squeeze(mdat(v,p,m,:))',...
%                 squeeze(mdat_ci(v,p,m,:,1))',squeeze(mdat_ci(v,p,m,:,2))');
            hp = plot(1:nbins,mdat,'.-');
            hp.MarkerSize = 20;
            hp.LineWidth = 2;
            hp.Color = cml(mi,:); mi = mi+1;
            
                
            clear mdat
        end
        
        % format & label axes
        ax = gca;
        ax.XTick = 1:nbins;
        if p == np            
            ax.XTickLabel = binlabels;
            ax.XTickLabelRotation = -65;
            ax.FontSize = 10;
        else
            ax.XTickLabel = [];
        end
        if v > 1
            ax.YTickLabel = [];
        end
        if p == 1
            title(voilist{v});
        end
        if v == 1
            ylabel([parlabels{p}]);
            if p == np                
                xlabel('distance from attention target');
            end
        end
        
        % do some last minute axis stuff
        pp = pp+1; xlim([0 7]); ylim(yl(p,:));
    end
end
hl = legend({'real data',mlabels{plotm}});
hl.Box = 'off';
hl.Position = [0.0197 0.8492 0.0650 0.0856]; % taken from manual adj