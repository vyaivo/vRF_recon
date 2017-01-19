function vRFRecon_extractTrialData(subList,voi_names)
% This function reads in the BOLD data files (already passed through
% compVOI to crop out voxels not in localizer). It compiles them into a
% data structure with all the relevant information needed to either
% calculate voxel receptive fields or do stimulus reconstructions in one of
% the two saved data structures,

% The mrdat struct contains all the MR timecourses and metadata. Everything
% is in TR multiples.
% The trialdat struct contains estimated beta weights for each trial, along
% with attention condition labels, stimulus info, and run labels. This
% should be enough to do vRF or recon analysis -- you'll just need to
% separate training fr. testing runs.

% VAV 6/14/2016
% updated VAV 12/8/2016

%% argument parsing
if nargin < 1
    subList = { {'AA91','AA92','AA93'}, {'AI92','AI93','AI94'}, ...
        {'AL91','AL92','AL93'}, {'AP91','AP92','AP93'}, ...
        {'AR91','AR92','AR93'},{'AT91','AT92','AT93'}, {'AU91','AU92','AU93'} };
    voi_names = {'V1','V2','V3','V3AB','V4','IPS0'};
elseif nargin < 2
    voi_names = {'V1','V2','V3','V3AB','V4','IPS0'};
end

%% path vars
root = '/usr/local/serenceslab/vy/vRFRecon_OSF/';
datadir = 'vRFRecon_mrStructs';
trialdir = 'vRFRecon_trialData';
if ~exist(fullfile(root,datadir),'dir')
    mkdir(fullfile(root,datadir));
end

runTypes = {'JitterLeft','JitterRight','JitterFix','DiscreteLeft',...
    'DiscreteRight'};

hemis = {'LH','RH'};
locstr = 'locSquareFull';

%%
for v = 1:length(voi_names)
    voi = voi_names{v};
    
    for s = 1:length(subList)  % loop over unique subjects
        sub = subList{s}{1}(1:2);        
        
        for tr = 1:length(runTypes)
%             for aa = 1:length(ttlist{tt})                
%                 fn2 = cell(2,length(subList{s}));
                for sx = 1:length(subList{s})
                    sn = subList{s}{sx};
                    for hh = 1:2        % combine data from both hemispheres
                        fn2{hh,sx} = sprintf('%s%s/%s_%s_Norm2_%s_%s-%s.mat',...
                            root,datadir,sn,locstr,runTypes{tr},hemis{hh},...
                           voi_names{v});
                    end
                end
                [mrdat.(runTypes{tr}), trialdat.(runTypes{tr})] = ...
                    getSessData(fn2);
%             end
        end

        savefn = sprintf('%s%s/%s_%s_vRFRecon_allTrialData.mat',root,...
            trialdir,sub,voi);
        save(savefn,'mrdat','trialdat','-v7.3');
        
        fprintf('SAVED %s\n',savefn);
        clear mrdat trialdat
    end
end

end     % end func

function [mrcat, b_resp] = getSessData(sesslist)
% subfunction to pull out data and combine across sessions
% The input variable 'sesslist' is a list of filenames that should be
% compiled into one data structure. This allows you to pass in one file
% list for a desired training set, one for a desired test set, etc...

% some options when loading the fMRI data from the BV structs
hemis = {'LH','RH'};

mrcat = struct(hemis{1},[],hemis{2},[],'all_tc',[],'attn_label',[],...
    'run_label',[],'prt_name',[]);
b_resp = struct('betaWeights',[],'stimStruct',[],'stimPos',[],'R2',[],...
    'goodVox',[],'attn_label',[],'run_label',[]);

lastrun = 0;
%%
for ss = 1:size(sesslist,2)         % loop across sessions
    
    tmptc = [];
    for hh = 1:2                    % combine data across hemispheres
        %% load data
        fn = sesslist{hh,ss};

        if ~exist(fn,'file')
            fprintf('MISSING file %s!\n', fn);
            break;
        end
        fprintf('Loading %s...\n', fn);
        load(fn);
        
        if ss > 1 && ~isempty(mrcat.(hemis{hh}))
            mrcat.(hemis{hh})(ss) = mr;
        else
            mrcat(1).(hemis{hh}) = mr;
        end
        
        % save out timecourse: nTRs x nvoxels (combining L/R hemis)
        tmptc = cat(2,tmptc,mr.voiData(1).mod);
        
        if hh == 1
            % these are the same across hemispheres, so just combine here
            % across session
            
            % save out condition label
            if ~isempty(strfind(fn,'Left'))
                alab = 1*ones(size(mr.voiData(1).mod,1),1);
            elseif ~isempty(strfind(fn,'Right'))
                alab = 2*ones(size(mr.voiData(1).mod,1),1);
            elseif ~isempty(strfind(fn,'Fix'))
                alab = 3*ones(size(mr.voiData(1).mod,1),1);
            end
            mrcat.attn_label = cat(1,mrcat.attn_label,alab);
            clear alab;
            
            % save out prt label
            mrcat.prt_name = cat(2,mrcat.prt_name,mr.prtName);
            
            % save out run label
            for rr = mr.r
                mrcat.run_label = cat(1,mrcat.run_label,repmat(mr.r(rr)+lastrun,...
                    size(mr.voiData(1).mod,1)/(max(mr.r)),1) );
            end
            lastrun = max(mr.r) + lastrun;
        end

        clear mr
    end     % end hemisphere loop    

    if isempty(mrcat.run_label)
        % file DNE (line 143); keep going
        continue;
    end
        
    % nTRs x nvoxels (combining L/R hemis)
    mrcat.all_tc = cat(1,mrcat.all_tc,tmptc);

end     % end session loop

%% now calculate trial-by-trial beta weights using a GLM

for ri = 1:length(mrcat.prt_name)
    [bEst,stimdat,stimPos,meanr2,voxGood] = extractTrialBetas(...
        mrcat.all_tc(mrcat.run_label==ri,:),mrcat.prt_name{ri},...
        mrcat.LH(1).NofTRs,mrcat.LH(1).TR);
    b_resp.betaWeights = cat(1,b_resp.betaWeights,bEst);
    b_resp.stimStruct = cat(1,b_resp.stimStruct,stimdat);
    b_resp.stimPos = cat(1,b_resp.stimPos,stimPos);
    b_resp.R2 = cat(1,b_resp.R2,meanr2);
    b_resp.goodVox = cat(1,b_resp.goodVox,voxGood);
    b_resp.attn_label = cat(1,b_resp.attn_label,repmat(...
        mrcat.attn_label(find(mrcat.run_label==ri,1)),size(bEst,1),1));
    b_resp.run_label = cat(1,b_resp.run_label,repmat(ri,size(bEst,1),1));
    clear bEst stimdat stimPos meanr2 voxGood
end

end         % end func