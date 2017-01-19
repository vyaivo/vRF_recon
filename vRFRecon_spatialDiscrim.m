function par_ds = vRFRecon_spatialDiscrim(parcombs,plab)

% Calculates spatial discriminability metric (similar to Fisher info) based
% on the sum of spatial derivatives of vRFs in a defined region of space.
% Normalizes by the maximum of the group of derivatives so that the
% discriminability values are comparable across spatial position despite
% changes in vRF coverage, vRF responsiveness, etc.
% VAV 9/29/2016
% cleaned up for OSF 12/13/2016

if nargin < 1
    parcombs = [1 1 1 1 1;
                1 1 0 1 0;
                1 1 1 0 0;            
                0 0 1 1 0;
                1 1 0 0 0;
                0 0 0 1 0;
                0 0 1 0 0;
                0 0 0 0 0];
end
if nargin < 2    
    plab = {'full','p+a','p+s','s+a','p','s','a','none'};
end

%%
root = load_root;
fdir = 'vRFfits';
vdir = 'vRFs';
rdir = 'recons';
dsdir = 'spatialDiscrim';
skipsub = [];

clist = {'Left','Right','Fix'};

% load vrf data
vrfdat = fullfile(root,fdir,'AllSubs_CompiledvRFs.mat');
load(vrfdat);
% get outlier vRFs
load(sprintf('%s%s/all_vRF_diffscores_noOutliers.mat',root,fdir),'allout');

% define some params based on the first sub/VOI
exampleFile = sprintf('%s%s/%s_TrnAllJitter_acrossSess_V1.mat',...
    root,rdir,sublist{1});
load(exampleFile,'chanX','chanY');

% define the resolution of the simulated vRF activity using first sub
% this should be much smaller than the IEM matrix
stepsize = 0.1;
xlist = min(chanX):stepsize:max(chanX);
ylist = min(chanY):stepsize:max(chanY);    
[rx, ry] = meshgrid(xlist,ylist);
xx = rx(:);
yy = ry(:);

savefn = fullfile(root,dsdir,'vRFSpatialDiscrim.mat');

%%
for voin = 1:length(voilist)
    for sn = 1:length(sublist)

        subdat = alldat{sn,voin};
        % remove outliers
        subdat(allout{sn,voin},:,:) = [];

        if isempty(subdat) || size(subdat,1) == 1 || ismember(sn,skipsub)
            % if there's not enough data or we should skip the subject, skip
            disp('not enough data');
            continue;
        end
        nvox = size(subdat,1);        

        %% loop through & calculate the discriminability metric

        for p = 1:size(parcombs,1)

            vrf_dat = subdat;
            nochange = find(~parcombs(p,:));
            for ni = 1:length(nochange)
                vrf_dat(:,1:2,nochange(ni)) = repmat(vrf_dat(:,3,nochange(ni)),1,2,1);
            end
            
            vrfw = nan(nc,nvox,numel(xx));
            
            % loop through each condition & voxel    
            for cc = 1:nc
                for vv = 1:nvox

                    % make vRF weight matrix
                    vrfw(cc,vv,:) = vrf_dat(vv,cc,5) + make2dcos(xx,yy,...
                        vrf_dat(vv,cc,1),vrf_dat(vv,cc,2),...
                        vrf_dat(vv,cc,3),7)*vrf_dat(vv,cc,4);
                    im1 = reshape(vrfw(cc,vv,:),numel(ylist),numel(xlist));

                    % take the spatial derivative of the vRF
                    for yi = 2:numel(ylist)
                        dxw(yi-1,:) = (diff(im1(yi,:))/stepsize).^2;
                    end
                    for xi = 2:numel(xlist)
                        dyw(:,xi-1) = (diff(im1(:,xi))/stepsize).^2;
                    end

                    % sum across x/y dimension
                    all_ds(cc,vv,:,:) = dxw + dyw;

                end
            end
            sum_ds = squeeze(sum(all_ds,2));

            for c = 1:3
                % normalize each measurement by the max derivative value
                normds(c,:,:) = sum_ds(c,:,:) ./ max(max(sum_ds(c,:,:)));
            end

            par_ds(sn,voin,p,:,:,:) = normds;
            
            clear im1 vrfw all_ds norm_ds dxw dyw
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% make plots
 
save(savefn,'par_ds','parcombs','sublist','voilist','xx','yy','plab');

%%
% sum across subs
subds = squeeze(nanmean(par_ds));
ns = length(sublist);
atlocs = cat(1,atlocs,[0.0000    0.0000]);
np = length(plab);

%% just look at area near attended location
niters = 10000;
rs = 3582339;
rng(rs);

for v = 1:nv
    disp(['Doing ',voilist{v}]);
    for p = 1:np
        tmp = []; tmp2 = []; sl = [];
        
        for c = 1:2
            ic = setxor(c,[1 2]);
            [xi,yi] = find( sqrt( (atlocs(c,1)-rx).^2 + (atlocs(c,2)-ry).^2 ) <= 1 );
            for ii = 1:length(xi)
                tmp = cat( 1,tmp,squeeze(par_ds(:,v,p,c,xi(ii),yi(ii))) );
                tmp2 = cat( 1,tmp2,squeeze(par_ds(:,v,p,ic,xi(ii),yi(ii))) );
                sl = cat( 2,sl,1:ns );
            end
            clear xi yi
        end
        
        discrim_at(v,p) = mean(tmp(:));
        discrim_ig(v,p) = mean(tmp2(:));
        ati = nan(1,niters); igi = ati;
        parfor i = 1:niters
            ats = nan(1,ns); igs = nan(1,ns);
            for s = 1:ns
                % shuffle within subs
                si = find(sl==s);
                ats(s) = mean( tmp( si(randi(numel(si),1,numel(si))) ) );
                igs(s) = mean( tmp2( si(randi(numel(si),1,numel(si))) ) );
            end
            ati(i) = mean(ats);
            igi(i) = mean(igs);
        end
        allat{v,p} = tmp;
        allig{v,p} = tmp2;
        ciat(v,p,:) = prctile(ati,[2.5 97.5]);
        ciig(v,p,:) = prctile(igi,[2.5 97.5]);
        pval_atMig(v,p) = 2*min([mean(ati-igi<0),mean(ati-igi>0)]);
        
    end
end

[fv,fdrmask] = fdr(pval_atMig,0.05);

%%
save(savefn,'par_ds','parcombs','sublist','voilist','xx','yy',...
    'plab','ns','nv','np','discrim_at','discrim_ig','ciat','ciig','niters',...
    'rs','atlocs','clist','pval_atMig','fdrmask','-v7.3');
disp('SAVED.');

%% plot avg discrim

cm = lines;
figure;
h1 = bar(1:2:nv*2,discrim_at(:,1)',0.4,'FaceColor',cm(1,:)); hold on;
h2 = bar(2:2:nv*2,discrim_ig(:,1)',0.4,'FaceColor',cm(2,:)); hold on;
h1.EdgeColor = 'none'; h2.EdgeColor = 'none';
h2(1).Parent.XTick = 1.5:2:nv*2;
h2(1).Parent.XTickLabel = voilist;
h3 = ploterr(1:2:nv*2,discrim_at(:,1)',[],{ciat(:,1,1),ciat(:,1,2)},'k.','hhy',0.5);
h4 = ploterr(2:2:nv*2,discrim_ig(:,1)',[],{ciig(:,1,1),ciig(:,1,2)},'k.','hhy',0.5);
h3(1).Marker = 'none'; h4(1).Marker = 'none';
h3(2).LineWidth = 1; h4(2).LineWidth = 1;
xlim([0 nv*2+1]);
ylabel('Fine spatial discriminability');
box off;
legend('attend','ignore');

%% just plot for attend condition
cm2 = lines(6);

figure;
plotc = np:-1:2;
for v = 1:nv
    hh2 = ploterr( (v-1)*(np-1)+1:1:(v-1)*(np-1)+np-1, discrim_at(v,plotc),...
        [], {ciat(v,plotc,1),ciat(v,plotc,2)});
    hh2(1).LineWidth = 2; hh2(1).Color = cm2(v,:);
    hh2(2).LineWidth = 2; hh2(2).Color = cm2(v,:);
    box off; hold on;
end
labs = repmat(plab(2:end),1,nv);
hh2(1).Parent.XTick = 1:(np-1)*nv;
hh2(1).Parent.XTickLabel = flipud(horzcat(labs(:)));
hh2(1).Parent.XTickLabelRotation = 40;
xlim([0 (np-1)*nv+1]);
legend(voilist);
ylabel('Fine spatial discriminability');