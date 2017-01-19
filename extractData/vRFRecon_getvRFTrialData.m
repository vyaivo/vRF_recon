% vRFRecon_getvRFTrialData.m
% VAV 12/8/2016

subList = {'AA','AI','AL','AP','AR','AT','AU'};
voi_names = {'V1','V2','V3','V3AB','V4','IPS0'};

%% path vars
root = load_root;
trialdir = 'vRFRecon_trialData';
if ~exist(fullfile(root,trialdir),'dir')
    mkdir(fullfile(root,trialdir));
end

% we are calculating vRFs separately for each condition, so separate them
% as such...
condRuns = {{'JitterLeft','DiscreteLeft'},{'JitterRight','DiscreteRight'},...
    {'JitterFix'}};
nconds = length(condRuns);
hemis = {'LH','RH'};
locstr = 'locSquareFull';

%%
for v = 1:length(voi_names)
    
    for s = 1:length(subList)  % loop over subj
        
        %%
        fn = sprintf('%s%s/%s_%s_vRFRecon_allTrialData.mat',root,trialdir,...
            subList{s},voi_names{v});
        load(fn,'trialdat','mrdat');
        
        % need to balance data evenly across L, R, and fix
        nruns = zeros(nconds,1);
        for cc = 1:nconds
            tmprunno = 0;
            for rt = 1:length(condRuns{cc})
                tmprunno = max(trialdat.(condRuns{cc}{rt}).run_label) + tmprunno;
                if ~exist('tperrun','var')
                    tperrun = numel(trialdat.(condRuns{cc}{rt}).run_label) / tmprunno;
                end
            end
            nruns(cc) = tmprunno;
        end
        nBalanceRuns = min(nruns);
        shufruns = nan(nconds,nBalanceRuns);
        
        %%
        rs = 658912 + s + v;
        rng(rs);
        
        for cc = 1:nconds
            
            tmpb = []; tmppos = []; tmpr2 = []; rlab = [];
            lastr = 0;
            for rt = 1:length(condRuns{cc})
                if trialdat.(condRuns{cc}{rt}).attn_label(1) == cc
                    tmpb = cat(1,tmpb,...
                        trialdat.(condRuns{cc}{rt}).betaWeights);
                    tmppos = cat(1,tmppos,...
                        trialdat.(condRuns{cc}{rt}).stimStruct);
                    tmpr2 = cat(1,tmpr2,trialdat.(condRuns{cc}{rt}).R2);
                    rlab = cat(1,rlab,trialdat.(condRuns{cc}{rt}).run_label+lastr);
                    lastr = max(rlab);
                end
            end
            
            bdat{cc} = tmpb;
            stimDat{cc} = tmppos;
            R2{cc} = tmpr2;
            runLabel{cc} = rlab;
            
            % runs to grab when balancing datasets
            shufruns(cc,:) = randperm(nruns(cc),nBalanceRuns);
        end
        
        %% to enable L/R hemisphere analysis, label voxels as L or R
        nvox = size(bdat{1},2);
        hemiLabel = [];
        for hh = 1:length(hemis)
            hemiLabel = cat(1,hemiLabel,repmat(hh,...
                mrdat.DiscreteLeft.(hemis{hh})(1).voiData.size,1));
        end
        if numel(hemiLabel) ~= nvox
            fprintf('ERROR: mismatching voxel hemisphere labels. RETURNING...\n');
            return;
        end
        
        %% save data
        savefn = sprintf('%s%s/%s_vRF_acrossSess_wBalanceInds_%s.mat',...
            root,trialdir,subList{s},voi_names{v});
        save(savefn,'bdat','stimDat','R2','runLabel','hemiLabel','nvox',...
            'shufruns','rs');
        
        fprintf('SAVED %s\n',savefn);
        
        clear bdat stimDat R2 runLabel shufruns rs
    end
end