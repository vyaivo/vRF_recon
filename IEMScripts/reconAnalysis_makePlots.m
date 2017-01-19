%%% reconAnalysis_makePlots
%%% VAV 4/23/2015
%%% cleaned up for OSF 12/14/2016

%% path stuff, and load the data
root = load_root;
fitdir = 'reconfits';
trnfix = 'TrnAllJitter';

reconfn = sprintf('%s%s/AllSubs_CompiledRecons_%s.mat',root,fitdir,trnfix);
load(reconfn,'allRecons','reconSubAvg','reconFits','attendlocs','locs',...
    'nv','ns','voilist','sublist','res','fov');

%%
% for recon plots, figure out where to plot each recon so that the subplot
% grid follows the stimulus position grid
newlocs(:,1) = (locs(:,1)+(fov(1)/2)) * (res(1)/fov(1));
newlocs(:,2) = (locs(:,2)+(fov(2)/2)) * (res(2)/fov(2));
poslocs(:,1) = (newlocs(:,1)/res(1)) - 0.0025;
poslocs(:,2) = 1-(newlocs(:,2)/res(2));

% Figure 4B - average recon for AI, V1 attend left
sn = 2; voin = 1; c = 1;
if exist('allRecons','var')
    % save mem
    subdat = squeeze(allRecons(sn,voin,:,:,:));
    clearvar allRecons;
end
figure;
for p = 1:size(subdat,2) 
    h(p) = subplot('Position',[poslocs(p,1) poslocs(p,2) 0.05 0.05]);
    imagesc(reshape(subdat(c,p,:),res(2),res(1)),[-0.25 1.25]);
    axis equal; axis off;
end
suptitle('AI V1: attend left');            
match_clim(h);
cb = colorbar;
cb.Position = [0.9 0.1 0.025 0.75];
clear h;

%% make plots of all subject averaged recons
% Figure S4 - all averaged recons

% sort the locations
locs = round(locs,1);
leftlocs = find(locs(:,1) < 0);
% left & right are not matched. match 'em
for i = 1:length(leftlocs)
    rightlocs(i) = find(locs(:,1) == -1*locs(leftlocs(i),1) & ...
        locs(:,2) == locs(leftlocs(i),2) );
end
catloc = [leftlocs'; rightlocs];

flipbool = 1;

for v = 1:length(voilist)
    if flipbool
        % now flip so that left is always attended
        for c = 1:2
            ic = setxor(c,1:2);
            if c == 1
                atd(c,:,:) = reconSubAvg(v,c,catloc(c,:),:);
                igd(c,:,:) = reconSubAvg(v,c,catloc(ic,:),:);
            else
                % flip all of these...
                for p = 1:size(catloc,2)
                    atd(c,p,:) = reshape( fliplr(reshape(...
                        squeeze(reconSubAvg(v,c,catloc(c,p),:)),res(2),res(1))), ...
                        1,res(2)*res(1) );
                    igd(c,p,:) = reshape( fliplr(reshape(...
                        squeeze(reconSubAvg(v,c,catloc(ic,p),:)),res(2),res(1))), ...
                        1,res(2)*res(1) );
                end
            end
        end
        fliprec = cat(2,atd,igd);
        flipavg = squeeze(mean(fliprec));

        figure;pp=1;        
        for c = 1:2
            for p = 1:size(catloc,2)
                h(p) = subplot('Position',[poslocs(catloc(c,p),1) ,...
                    poslocs(catloc(c,p),2), 0.05, 0.05]);
                imagesc(reshape(flipavg(pp,:),res(2),res(1)));
                axis equal; axis off;
                pp = pp+1;
            end
        end
        suptitle(voilist{v});
        match_clim(h);
        cb = colorbar;
        cb.Position = [0.9 0.1 0.025 0.75];

        clear h cb atd igd fliprec flipavg
    else
        clist = {'left','right'};
        for c = 1:2
            figure;
            for p = 1:size(reconSubAvg,3) 
                h(p) = subplot('Position',[poslocs(p,1) poslocs(p,2) 0.05 0.05]);
                imagesc(reshape(reconSubAvg(v,c,p,:),res(2),res(1)));
                axis equal; axis off;
            end
            suptitle([voilist{v} ': attend ' clist{c}]);            
            match_clim(h);
        end        
        cb = colorbar;
        cb.Position = [0.9 0.1 0.025 0.75];
    end
end

%% make plot of all positions w/ amp as intensity

% Figure 4C - size/amp plot
ddim = size(reconFits);

pmap = parula(1000);
arange = linspace(0.5,1.25,1000);       % amplitude range
srange = linspace(2,3.5,300);           % min / max size limits
% sort the locations
locs = round(locs,1);
leftlocs = find(locs(:,1) < 0);
% left & right are not matched. match 'em
for i = 1:length(leftlocs)
    rightlocs(i) = find(locs(:,1) == -1*locs(leftlocs(i),1) & ...
        locs(:,2) == locs(leftlocs(i),2) );
end

for v = 1:nv
    
    subplot(2,ceil(nv/2),v);
    axis equal;    
    ylim([-3 3]); xlim([-6 6]);
    if v == 1
        set(gca,'Units','Points');
        axpos = get(gca,'Position');
        set(gca,'Units','Normalized');
        szratio = axpos(3)/diff(xlim);
    end
    
    % pull out amplitudes at each position
    allamps = squeeze(mean(reconFits(v,:,:,:,4),2));
    % attended amplitudes
    amps(1,:) = mean(cat(1,allamps(1,leftlocs),allamps(2,rightlocs)));
    amps(2,:) = mean(cat(1,allamps(2,leftlocs),allamps(1,rightlocs)));
    % now find which index is closest to the amp value
    [~,~,ampdiv] = histcounts(amps,arange);
    
    % pull out sz at each position
    allsz = squeeze(mean(rad2fwhm(reconFits(v,:,:,:,3)),2));
    % attended sz
    szs(1,:) = mean(cat(1,allsz(1,leftlocs),allsz(2,rightlocs)));
    szs(2,:) = mean(cat(1,allsz(2,leftlocs),allsz(1,rightlocs)));
    % now find which index is closest to the sz value
    [~,~,szdiv] = histcounts(szs,srange);
    szarea = pi*(szs.^2);
    
    plot(0,0,'k*'); hold all;
    scatter(locs(leftlocs,1),-1*locs(leftlocs,2),szratio*szarea(1,:),pmap(ampdiv(1,:),:),'filled');
    scatter(locs(rightlocs,1),-1*locs(rightlocs,2),szratio*szarea(2,:),pmap(ampdiv(2,:),:),'filled');    
    
%     scatter(locs(leftlocs,1),-1*locs(leftlocs,2),szdiv(1,:),pmap(ampdiv(1,:),:),'filled');
%     scatter(locs(rightlocs,1),-1*locs(rightlocs,2),szdiv(2,:),pmap(ampdiv(2,:),:),'filled');    
    plot(attendlocs(:,1),-1*attendlocs(:,2),'rx');
    
    title(voilist{v});
    axis equal;    
    ylim([-3 3]); xlim([-6 6]);
    
    if v == nv
        caxis([0.5 1.25]);
        b = colorbar;
        b.Label.String = 'Amplitude (BOLD z-score)';        
    end
end
set(b,'Position',[0.92 0.1 0.03 0.8]);
suptitle('Averaged reconstruction size & amplitude at all positions, left is attended');

%% make a vector plot

% make a center shift vector plot
thstep = pi/6;
tharray = -pi:thstep:(pi-thstep);
ntbins = length(tharray)-1;
rstep = 0.8;
rarray = 0:rstep:4-rstep;
nrbins = length(rarray)-1;

[t r] = meshgrid(tharray,rarray);
% we have to rebin the data.
for v = 1:nv
    for c = 1:2
        ic = setxor(c,[1 2]);
        
        % convert x & y to theta & rho (centered to dist from attn loci)
        % now we have attend
        theta(:,1,:) = atan2( reconFits(v,:,c,:,2)-attendlocs(c,2), ...
            reconFits(v,:,1,:,1)-attendlocs(c,1) );  % theta = atan2(y,x);
        rho(:,1,:) = sqrt( (reconFits(v,:,c,:,1)-attendlocs(c,1)).^2 + ...
            (reconFits(v,:,1,:,2)-attendlocs(c,2)).^2 );
        % now let's do ignoresh
        theta(:,2,:) = atan2( reconFits(v,:,c,:,2)-attendlocs(ic,2), ...
            reconFits(v,:,1,:,1)-attendlocs(ic,1) );
        rho(:,2,:) = sqrt( (reconFits(v,:,c,:,1)-attendlocs(ic,1)).^2 + ...
            (reconFits(v,:,1,:,2)-attendlocs(ic,2)).^2 );
        
        tr_attend = [reshape(theta(:,1,:),ns*51,1) reshape(rho(:,1,:),ns*51,1)];
        tr_ignore = [reshape(theta(:,2,:),ns*51,1) reshape(rho(:,2,:),ns*51,1)];
        
        % now let's bin this!
%         [na,ca] = hist3(tr_attend,[tharray; rarray]);
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
                end
                clear inds
            end
        end

    end
end

allorigs = [(r(:).*cos(t(:)))+atlocs(c,1) (r(:).*sin(t(:)))+atlocs(c,2)];
% save out these points
binpts(c,:,:) = allorigs;
if plotbool
    subplot(1,2,c);
    h = quiver(allorigs(:,1),allorigs(:,2),xcomps',ycomps',...
        'LineWidth',2,'MaxHeadSize',0.25,'AutoScale','on',...
        'AutoScaleFactor',1.25);
    hold all;
    plot(0,0,'k.');
    plot(atlocs(c,1),atlocs(c,2),'ro');
    plot(atlocs(ic,1),atlocs(ic,2),'o','MarkerEdgeColor',[0.4 0.4 0.4]);
    axis equal;
    xlim([-6 6]); ylim([-4 4]);
end