% vRFRecon_getReconTrialData.m
% VAV 12/8/2016

subList = {'AA','AI','AL','AP','AR','AT','AU'};
voi_names = {'V1','V2','V3','V3AB','V4','IPS0'};

%% path vars
root = load_root;
trialdir = 'vRFRecon_trialData';
if ~exist(fullfile(root,trialdir),'dir')
    mkdir(fullfile(root,trialdir));
end

trnRuns = {'JitterLeft','JitterRight','JitterFix'};	% run names for training set
tstRuns = {'DiscreteLeft','DiscreteRight'};         % run names for testing set
ttlist = {tstRuns,trnRuns};                         % entire list of run names

hemis = {'LH','RH'};
locstr = 'locSquareFull';

%%
for v = 1:length(voi_names)
    
    for s = 1:length(subList)  % loop over subj
        
        rs = 46789252 + s + v;
        rng(rs);
        
        %%
        fn = sprintf('%s%s/%s_%s_vRFRecon_allTrialData.mat',root,trialdir,...
            subList{s},voi_names{v});
        load(fn,'trialdat');
        
        %% training data
        
        % balance the training set across L/R runs (fixation always
        % included)
        for tt = 1:length(trnRuns)
            maxruns(tt) = max(trialdat.(trnRuns{tt}).run_label);
            if ~exist('tperrun','var')
                tperrun = numel(trialdat.(trnRuns{tt}).run_label) / maxruns(tt);
            end
        end
        nBalanceRuns = maxruns(3);      % include as many fix runs as you want, they're neutral
        
        trn_c_all = []; trnb = []; trn_stim_all = [];
        for tt = 1:length(trnRuns)
            
            if maxruns(tt) > nBalanceRuns   % need to discard some runs
                useruns = randperm(maxruns(tt),nBalanceRuns);
            else
                useruns = 1:length(trialdat.(trnRuns{tt}).stimStruct);
            end
            
            % save out the relevant data. conditions can all be
            % concatenated since we're training across condition
            for rr = useruns
                starti = (rr-1)*tperrun + 1;
                endi = rr*tperrun;
                trnb = cat(1,trnb,trialdat.(trnRuns{tt}).betaWeights(starti:endi,:));
                trn_c_all = cat(1,trn_c_all,...
                    trialdat.(trnRuns{tt}).attn_label(starti:endi));
                trn_stim_all = cat(1,trn_stim_all,trialdat.(trnRuns{tt}).stimStruct(rr));                
            end            
            
        end
        
        %% testing data
        
        tst_b_all = cell(2,1); tst_c_all = cell(2,1); tst_pos_all = cell(2,1);
        for tt = 1:length(tstRuns)
            tst_b_all{tt} = trialdat.(tstRuns{tt}).betaWeights;
            tst_c_all{tt} = trialdat.(tstRuns{tt}).attn_label;
            tst_pos_all{tt} = trialdat.(tstRuns{tt}).stimPos;
        end
        
        %% save data
        
        savefn = sprintf('%s%s/%s_recon_TrnAllJitter_acrossSess_Bilat-%s.mat',...
            root,trialdir,subList{s},voi_names{v});
        save(savefn,'trn_c_all','trnb','trn_stim_all',...
            'tst_b_all','tst_c_all','tst_pos_all','rs');
        fprintf('SAVED %s\n',savefn);
    end
end

%% concatenate into big ROIs
bigROIs = { {'V1','V2','V3','V4'}, {'V3AB','IPS0'} };
nc = 2;
bigROIname = {'superVis','superParietal'};

for bb = 1:length(bigROIs)     
    
    for s = 1:length(subList)
        savebfn = sprintf('%s%s/%s_recon_TrnAllJitter_acrossSess_Bilat-%s.mat',...
            root,trialdir,subList{s},bigROIname{bb});
        
        bigtrnb = []; big_tstb = cell(1,nc);
        for vb = 1:length(bigROIs{bb})
            loadfn = sprintf('%s%s/%s_recon_TrnAllJitter_acrossSess_Bilat-%s.mat',...
                root,trialdir,subList{s},bigROIs{bb}{vb});
            load(loadfn);
            % concatenate data
            bigtrnb = cat(2,bigtrnb,trnb);
            for c = 1:nc
                big_tstb{c} = cat(2,big_tstb{c},tst_b_all{c});
            end
        end
        trnb = bigtrnb;
        tst_b_all = big_tstb;
        
        save(savebfn,'trn_c_all','trnb','trn_stim_all','tst_b_all',...
            'tst_c_all','tst_pos_all');
        fprintf('Saved %s\n',savebfn);
        
        clearvars trnb tst_b_all trn_c_all trn_stim_all tst_c_all tst_pos_all
    end
end