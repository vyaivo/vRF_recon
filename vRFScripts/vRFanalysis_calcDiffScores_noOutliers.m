%%% vRFanalysis_calcDiffScores_noOutliers.m
%%% Loads compiled vRF fit parameters (across all subjects & ROIs) and
%%% calculates difference scores that describe voxel attentional
%%% modulations (attend target - attend fixation). 

%%% VAV 1/13/2017 - separated diff score calculation from stats script.
%%% VAV 1/16/2017 - added outrows param to make sure all data associated w
%%% an outlier voxel is removed.
%%% VAV 1/18/2017 - added 9th row w voxel eccentricity to fitmat

% load the vRF fit parameter data
root = load_root;
fitdir = 'vRFfits';
fn = fullfile(root,fitdir,'AllSubs_CompiledvRFs.mat');
load(fn);

parlabel = {'position','size','amplitude','baseline'};
savefn = sprintf('%s%s/all_vRF_diffscores_noOutliers.mat',root,fitdir);

%% sort the vRFs

% This makes a giant matrix ('fitmat') with vRF difference scores as the first column,
% and all subsequent columns are labels that describe properties of the
% vRF. Here is the content within each column with example possible values
% within each column in square brackets:

% col 1: data (e.g. all difference scores) [-0.1, 0.02, 0.2,...]
% col 2: subject [s1, s2, s3...]
% col 3: region of interest [V1, V2, V3...]
% col 4: bin number (by vRF's distance from the attended target) [1, 5, 3 ...]
% col 5: distance from the attended target (continuous val) [0.1, 1.1, 0.6...]
% col 6: parameter label (position, size, amplitude, baseline) [1, 4, 2...]
% col 7: contralateral vs. ipsilateral to attention target [1, 1, 2,...]
% col 8: left or right hemisphere [2, 2, 1...]
% col 9: eccentricity, or dist from fixation (continuous val) [0.4, 0.8, 0.2...]

%% 

% set up dist from attend bins (for col 4)
% voxels > 2.5 degs from attention target are either off the screen or
% start crossing into the other vis hemifield (and their attn modulations
% would be contaminated by other factors)
distbins = 0:0.25:2.5;

% create empty label vectors to fill in for each column of big matrix
adat = [];slab = [];vlab = [];dlab = [];
blab = [];clab = [];voxlab =[];hlab=[]; eccdat = [];

for s = 1:ns
    for v = 1:nv
        % read in data from each sub/ROI
        if isempty(alldat{s,v})
            % no voxels exist in this sub/ROI pair
            continue;
        end
        %%
        % thisdat is nvox x nconds x fit params, where
        % cond 1 = L, cond 2 = R, cond 3 = fix
        % and params are: 1 = x pos, 2 = y pos, 3 = size, 4 = amp, 5 =
        % baseline
        thisdat = alldat{s,v}(:,:,1:npar);
        nvox = size(thisdat,1);
        
        %% calculate attend - fix difference scores
        
        % prealloc array for all difference scores for each voxel, each
        % condition, and each parameter (collapsing x/y into a 'position'
        % parameter)
        diffscores = nan(nvox,2,size(thisdat,3)-1);
        
        % prealloc array for voxel labels as contra or ipsi to attended
        % target. note that each voxel will appear twice -- for example,
        % once for contra attend L, once for ipsi attend R.        
        contra_ipsi = nan(2,nvox);
        
        % prealloc array for the distance covariate (col 5 in giant fitmat
        % matrix)
        dist_at_fix = nan(nc-1,nvox);
        eccen = sqrt( thisdat(:,3,1).^2 + thisdat(:,3,2).^2 );
        
        % iterate over attend L/R
        for c = 1:nc-1
            
            % distance from attn target during fix cond
            dist_at_fix(c,:) = sqrt( (thisdat(:,3,1)-atlocs(c,1)).^2 + ...
                (thisdat(:,3,2)-atlocs(c,2)).^2 );
            
            % sort vRFs by distance from attention locus during fix cond
            [~,bi(c,:)] = histc(dist_at_fix(c,:),distbins);
            
            % distance from attn target during attend L/R cond
            dist_attend = sqrt( (thisdat(:,c,1)-atlocs(c,1)).^2 + ...
                (thisdat(:,c,2)-atlocs(c,2)).^2 );
            
            % now calc diff between position measures above to get an
            % overall score for how much closer the vRF moves to the attn
            % target
            diffscores(:,c,1) = dist_attend' - dist_at_fix(c,:);
            
            for p = 3:npar
                % calculate difference from fixation condition for all
                % other parameters (e.g. param 3 is size, param 4 is amp,
                % param 5 is baseline)
                diffscores(:,c,p-1) = thisdat(:,c,p) - thisdat(:,3,p);
            end
            
            % LABEL CONTRA VS IPSI
            % A) first get voxel hemisphere labels in all_hemilabs
            % B) then get attn labels. if they ARE OPPOSITE (e.g., hemi = 1
            % & attn = 2), you get CONTRA (1). if they ARE THE SAME you get
            % IPSI (2).
            ipsi = find(all_hemilabs{s,v} == c);
            contra = setxor(ipsi,1:nvox);
            
            contra_ipsi(c,contra) = 1;
            contra_ipsi(c,ipsi) = 2;
            
            clear dist_attend ipsi contra
        end
        
        % also save out columns of voxel number labels so we can use them
        % later to figure out which ones are outliers. we are making this
        % matrix the same size as the diffscores matrix.
%         voxnum = repmat([1:size(diffscores,1)]',1,2,npar-1);
        catvoxdat = cat(1,repmat([1:size(diffscores,1)]',1,npar-1),...
            repmat([1:size(diffscores,1)]',1,npar-1));
        
        % save out a matrix of the same size, but with LH/RH labels for
        % each voxel
        cathemi = cat(1,repmat(all_hemilabs{s,v},1,npar-1),...
            repmat(all_hemilabs{s,v},1,npar-1));
        
        % also adding an absolute eccentricity column
        catecc = cat(1,repmat(eccen,1,npar-1),repmat(eccen,1,npar-1));
        
        %% concatenate all the data
        
        % we want to fold the data together across conditions (attend L vs R)
        % so we should get 2x the vRFs
        catdiff = squeeze(cat(1,diffscores(:,1,:),diffscores(:,2,:)));
%         cathemi = squeeze(cat(1,hemidat(:,1,:),hemidat(:,2,:)));
%         catvoxdat = squeeze(cat(1,voxnum(:,1,:),voxnum(:,2,:)));
        catbin = [bi(1,:) bi(2,:)];
        catfdist = [dist_at_fix(1,:), dist_at_fix(2,:)];
        catcipsi = [contra_ipsi(1,:),contra_ipsi(2,:)];
        
        % all the associated labels for these voxels need to be the same
        % size as the diff score matrix after we collapsed across attention
        % conditions. e.g. need subject labels for all voxels in that
        % matrix, ROI labels, etc.
        sl = repmat(s,size(catdiff));                % subject label
        vl = repmat(v,size(catdiff));                % ROI label
        cl = repmat(catcipsi',1,size(catdiff,2));    % contra/ipsi label
        bl = repmat(catbin',1,size(catdiff,2));      % bin label
        dl = repmat(catfdist',1,size(catdiff,2));    % distance fr attend
        
        % save out all the labels for the ANOVA in giant columns
        slab = cat(1,slab,sl);
        vlab = cat(1,vlab,vl);
        blab = cat(1,blab,bl);
        dlab = cat(1,dlab,dl);
        clab = cat(1,clab,cl);
        hlab = cat(1,hlab,cathemi);
        adat = cat(1,adat,catdiff);
        voxlab = cat(1,voxlab,catvoxdat);
        eccdat = cat(1,eccdat,catecc);
        
        clearvars dist_at_fix voxnum hemidat bi
        clearvars catdiff cathemi catvoxdat catbin catfdist catcipsi catecc
        clearvars thisdat diffscores sl vl bl dl cl
    end
end

%% now stick these labels all together in big matrix
% update number of params, now that we've collapsed the x/y fit parameters
% into an overarching 'position' (or distance fr. attn targ) metric
npar = npar-1;
% also set the number of attn hemifield labels (contra vs ipsi)
na = numel(unique(clab(:)));

% set up parameter labels (1 - 4)
plab = nan(size(adat));
for p = 1:npar
    plab(:,p) = p;
end

% create the giant data matrix!
fitmat = [adat(:) slab(:) vlab(:) blab(:) dlab(:) plab(:) clab(:) hlab(:) eccdat(:)];
% fitmat is: data, sub label, VOI label, bin label, distance covariate,
% param label, contra/ipsi (attn side) label, hemisphere label,
% eccentricity

clearvars adat slab vlab blab dlab plab clab hlab eccdat
clearvars alldat all_hemilabs

% correspondingly, make a column of vox #s that has the same number of rows
% as fitmat. we'll need this for the next part, where we figure out which
% voxels are outliers.
voxlab2 = voxlab(:);

%% Get outlier indices if vox changes > 3SD from mean
% The outlier thresholding is based on the difference scores across
% the whole population (e.g., all conditions, all subjects, both
% hemispheres). So we need to first find the outliers in this dist, then save out
% identifiers which specify which subject/region of interest those outliers
% belong to. Those indices are saved out in the final analysis data file so
% they can be used to constrain subsequent analyses.

% prealloc a cell array of all outliers for each VOI
outi = cell(nv,1);
% prealloc a cell array that separates all voxels by subjects/VOI
allout = cell(ns,nv);
% prealloc an array that has all fitmat rows to remove
outrows = [];

for v = 1:nv
    
    % (1) IDENTIFY OUTLIERS. this has to be done separately for each type
    % of diff score parameter (e.g., pos changes are on a different scale
    % than baseline changes) and for each VOI (population of pos changes is
    % different in V1 than IPS0).
    outi{v} = [];
    for p = 1:npar
        % grab all the data for this VOI / fit param
        rowi = find(fitmat(:,3)==v & fitmat(:,6)==p);   % row indices from fitmat
        thisdat = fitmat(rowi,1);                       % actual diff scores
        % find standard deviation of this population
        sd = std(thisdat);
        sout = find(abs(thisdat - mean(thisdat)) > 3*sd);
        % save fitmat row indices for these outliers
        outi{v} = cat(1,outi{v},rowi(sout));
        
        clear rowi thisdat sd sout
    end
    
    % (2) SAVE OUT OUTLIERS IN allout{subj,voi}
    
    % figure out the subject labels for each of the outliers in voilabel == v
    % find unique sub labels for these outlier voxels
    thissubs = unique(fitmat(outi{v},2));
    % now for each subject, need to get the original voxel # label (e.g.,
    % vox 35 for sub 1, VOI 1; 
    % vox 9 for sub 2, VOI 1, etc.)
    for si = 1:length(thissubs)
        % grab fitmat row indices for this sub/VOI
        thiss = find(fitmat(:,2)==thissubs(si) & fitmat(:,3)==v);

        % now see if any of them match the voxel outliers we found (i.e. a
        % set intersection of all voxels for this sub/VOI and for all
        % outlier voxels for this VOI)
        subouti = intersect(outi{v},thiss);
        
        % these indices don't capture all of the rows of fitmat that fit
        % that voxel (just the ones that are for that param label.) So need
        % to get voxel labels and get those rows also
        thesevox = voxlab2(subouti);
        voxrowi = find(ismember(voxlab2,thesevox));
        
        % now get the specific voxel #s for this sub/VOI and save them
        allout{thissubs(si),v} = voxlab2(subouti);
        
        % also save out all the fitmat rows for this sub/VOI
        outrows = cat(1,outrows,voxrowi);
        
        clear thiss subouti thesevox voxrowi
    end
    
    clear thissubs
end
clear voxlab
voxlabs = voxlab2;
clear voxlab2;

%% remove the outliers for these analyses
newfmat = fitmat;
% allouti = vertcat(outi{:});     % this is all the outliers in the giant matrix fitmat
% newfmat(allouti,:) = [];        % take them out!
newfmat(outrows,:) = [];
%%
fitmat = newfmat;
clear newfmat allouti

%% SAVE EVERYTHING
clear s v c p si;   % clear the iterators from the workspace

save(savefn);
fprintf('SAVED %s\n',savefn);

%% clean up
clear;