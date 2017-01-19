% simAnalysis_calculateRMSEwFullData_bigROIs.m
% doing RMSE comparing with full dataset recons (e.g., all voxels, not just
% voxels w good vRFs). This script will also fit the results of the layered
% IEM / simulation if the doFit flag is set to 1. Finally, this script will
% save out the difference scores for the full data recons, the reduced data
% recons, and the modeled recons before calculating the RMSE.

% cleaned up for OSF 12/16/2016
% VAV 1/6/2017 - changed to accommodate big super ROIs

%%
% fit options
doFit = 0;      % fit the layered IEM recons
niters = 100;   % resample the data & fit for this many iterations

% get some data in place
root = load_root;
simdir = 'sims';

% load full recon fits & some other accompanying variables
fitdir = 'reconfits';
fitfix = 'acrossSess';
trnFix = 'goodvRFs_TrnAllJitter';
fn = sprintf('%s%s/AllSubs_CompiledRecons_bigROIs_%s.mat',root,fitdir,trnFix);
load(fn,'sublist','locs','attendlocs','fov','res','voilist','nv','nbins','nc','npar');
atlocs = attendlocs;
np = npar - 1;
npos = size(locs,1);
clearvars attendlocs npar

% need this for calculating difference scores in recon fits later on
locs = round(locs,1);
condlocs(1,:) = find(locs(:,1) < 0);    % any x pos < 0 is attend left
for i = 1:size(condlocs,2)
    condlocs(2,i) = find(locs(:,1) == -1*locs(condlocs(1,i),1) & ...
        locs(:,2) == locs(condlocs(1,i),2) );
end

% name of giant file with all fit data
bigfitfn = sprintf('%s%s/all_sim_bigROIs_fitsWModel.mat',root,simdir);
% labels of each recon parameter (position,size,amplitude,baseline)
parlab = {'p','s','a','b'};

%% find out which subjects we had to exclude
fn = dir(sprintf('%s%s/sim_vRFRecon_allSubsNoOutliers_super*.mat',root,simdir));

for ff = 1:length(fn)
    load(fn(ff).name,'skipsub','voin');
    badsub = unique(skipsub);
    fprintf('ROI %d: skipped ', voin);
    for bb = 1:length(badsub)
        fprintf('%s ', sublist{badsub(bb)});
    end
    fprintf('\n');
    clear skipsub voin badsub
end

%% get fits & calculate difference scores
% 
if doFit
    simAnalysis_fitSimRecons(niters,voilist,fov);
else    
    if ~exist('bigfitfn','file')
        for v = 1:nv
            fitfn = sprintf('%s%s/sim_RMSE_fit_%s.mat',root,simdir,voilist{v});
            load(fitfn,'allfits','niters');

            if ~exist('parcombs','var')
                thisfn = sprintf('%s%s/sim_vRFRecon_allSubsNoOutliers_%s.mat',...
                    root,simdir,voilist{v});
                load(thisfn,'parcombs');
                nm = size(parcombs,1);
                mlabels = cell(1,nm);
                for mi = 1:nm
                    mlabels{mi} = horzcat(parlab{find(parcombs(mi,2:end))});
                end
                modeldat = nan(nm,nv,niters,nc,24,np);
            end

            tmp = reshape(allfits,niters,nm,nc,npos,[]);
            % now tmp is 100 x 11 x 2 x 51 x 5
            % convert size constant to FWHM
            tmp(:,:,:,:,3) = rad2fwhm(tmp(:,:,:,:,3));
            % compute difference scores for model fits
            for mi = 1:nm
                modeldat(mi,v,:,:,:,:) = recon_attendMinusIgnore_acrossBlocks(...
                    squeeze(tmp(:,mi,:,:,:,:)),atlocs,condlocs);
            end
            clear allfits tmp
        end

        % save into big matrix
        save(bigfitfn,'modeldat','parcombs','mlabels','niters','nm',...
            'voilist','nv','niters','nm','condlocs','nbins','locs',...
            'atlocs','-v7.3');
    else
        load(bigfitfn);
    end
end

%% load other data
fn2 = fullfile(root,'reconfits','AllSubs_CompiledRecons_bigROIs_TrnAllJitter.mat');
load(fn2,'reconFits');
reconFits(:,:,:,:,3) = rad2fwhm(reconFits(:,:,:,:,3));
reconFits = squeeze(nanmean(reconFits,2)); % avg across subs
fullrec = recon_attendMinusIgnore_acrossBlocks(reconFits,atlocs,condlocs);
clear reconFits;

% load real (reduced) recon fits
load(fn, 'reconFits');
reconFits(:,:,:,:,3) = rad2fwhm(reconFits(:,:,:,:,3));
reconFits = squeeze(nanmean(reconFits,2)); % avg across subs
reducedrec = recon_attendMinusIgnore_acrossBlocks(reconFits,atlocs,condlocs);
clear reconFits;

save(bigfitfn,'modeldat','fullrec','reducedrec','parcombs','mlabels',...
    'voilist','nv','niters','nm','condlocs','nbins','locs','atlocs','-v7.3');

%% calculate RMSE!
rs = 39502385;
rng(rs);

reducedRMSE = nan(nv,1);
ci_reducedRMSE = nan(nv,2);
bootRMSE = nan(nv,niters);
modelRMSE = nan(nv,nm);
ciRMSE = nan(nv,nm,2);

for v = 1:nv
    rdat = reshape(fullrec(v,:,:,:)-reducedrec(v,:,:,:),...
        size(fullrec,2)*size(fullrec,3),[]);
    for i = 1:niters
        ri = randsample(size(rdat,1),size(rdat,1),1);
        bootRMSE(v,i) = sqrt( nanmean( reshape(rdat(ri,:),1,[]).^2 ));
    end
    reducedRMSE(v) = nanmean(bootRMSE(v,:),2);
    ci_reducedRMSE(v,:) = abs(reducedRMSE(v) - prctile(bootRMSE(v,:),...
        [2.5,97.5],2));
    
    % now for the layered IEM RMSE
    for mi = 1:nm
        iterrmse = repmat(fullrec(v,:,:,:),niters,1,1,1) - ...
            squeeze(modeldat(mi,v,:,:,:,:));
        mi1 = sqrt( nanmean(nanmean(mean(iterrmse.^2,2),3),4) );
        modelRMSE(v,mi) = nanmean(mi1);
        ciRMSE(v,mi,:) = abs(modelRMSE(v,mi) - prctile(mi1,[2.5,97.5],1));
        clear iterrmse mi1
    end
    
end

% %%
% figure;
% h = errorbar(1:nv,reducedRMSE,ci_reducedRMSE(:,1),ci_reducedRMSE(:,2),'.-',...
%     'LineWidth',2);
% h.Parent.XTickLabel = voilist;
% 
% %%
% figure;
% for v = 1:nv
%     sh(v) = subplot(2,nv/2,v);
%     errorbar(1:nm,modelRMSE(v,:),ciRMSE(v,:,1),ciRMSE(v,:,2),'.-',...
%         'LineWidth',2); hold on;
%     plot(0,reducedRMSE(v),'ro','LineWidth',2); hold on;
%     sh(v).XTick = 0:nm;
%     sh(v).XTickLabels = ['real(reduced)', mlabels];
%     sh(v).XTickLabelRotation = 35;
%     box off;
%     title(voilist{v});
% end
% %match_ylim(sh); 
% match_xlim(sh,[-1 nm+1]);

%%
[~,bestm] = min(modelRMSE,[],2);
[~,worstm] = max(modelRMSE,[],2);
% now test if they're significantly different
for v = 1:nv
    iterrmse = repmat(fullrec(v,:,:,:),niters,1,1,1) - ...
            squeeze(modeldat(nm,v,:,:,:,:));
    mi2 = sqrt( mean(mean(mean(iterrmse.^2,2),3),4) ); clear iterrmse;
    for m = 1:nm-1
        iterrmse = repmat(fullrec(v,:,:,:),niters,1,1,1) - ...
            squeeze(modeldat(m,v,:,:,:,:));
        mi1 = sqrt( mean(mean(mean(iterrmse.^2,2),3),4) );        
       
        pval_eachm(v,m) = 2*min([mean(mi1-mi2 < 0),...
                mean(mi1-mi2 > 0)]);
    end
end

[f,fmask] = fdr(pval_eachm,0.05);

%%
savefn = sprintf('%s%s/layeredSimModel_bigROIs_fitRMSE.mat',root,simdir);
save(savefn,'modelRMSE','ciRMSE','rs','reducedRMSE','voilist','nv',...
    'ci_reducedRMSE','bootRMSE','pval_eachm','fmask','parcombs','mlabels');

%% Figure 6C

sorti = [11 8 7 5 3 2 1 10 6 4 9];

figure;
fm = fmask(:,sorti(2:end));
peach = pval_eachm(:,sorti(2:end));

for v = 1:nv
    sh(v) = subplot(2,nv/2,v);
    errorbar(1:nm,modelRMSE(v,sorti)-modelRMSE(v,nm),...
        ciRMSE(v,sorti,1),ciRMSE(v,sorti,2),'.-',...
        'LineWidth',2); hold on;
    hh = refline(0,0);
    hh.Color = 'k';
    for m = 1:nm-1
        if fm(v,m)==1
            plot(m+1,0.05 ,'k*');
        elseif peach(v,m)<.05
            plot(m+1,0.05,'k+');
        end
    end
    sh(v).XTick = 1:nm;
    sh(v).XTickLabels = mlabels(sorti);
    sh(v).XTickLabelRotation = 35;
    box off;
    title(voilist{v});
end
match_ylim(sh); % -0.0500    0.0600
match_xlim(sh,[0 nm+1]);

%% Make Table 4 - mean + 95 CIs on RMSE
mlabels{end} = 'none';
modelRMSE = round(modelRMSE,3);
ci = ciRMSE;
ci(:,:,1) = modelRMSE - ci(:,:,1);
ci(:,:,2) = modelRMSE + ci(:,:,2);
ci = round(ci,3);
x = arrayfun(@cat,repmat(2,nv,nm),modelRMSE,ci(:,:,1),ci(:,:,2),...
    'UniformOutput',0);
tab4 = cell2table(x,'VariableNames',mlabels);
tab4.realData = cat(2,round(reducedRMSE,3),...
    round(reducedRMSE-ci_reducedRMSE(:,1),3),...
    round(reducedRMSE+ci_reducedRMSE(:,2),3));