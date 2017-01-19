function [bEst,stim,stimPos,mean_r2,voxGood] = extractTrialBetas(data,...
    pname,nTRs,TR)
% Extract the signal in each voxel on each trial, using trial-by-trial GLM.
% also calculates how well the GLM predicts the real data (R^2)
% updated/simplified/commented VAV 12/8/2016

% Define some HRF params
params.hpttp  = 5;      %     HRF time to peak for positive response (5)
params.hnttp = 15;      %     HRF5 time to peak for negative response (15)
params.hpnr = 6;        %     HRF positive/negative ratio (6)
params.hons = 0;        %     HRF onset (0)
params.hpdsp = 1;       %     HRF positive response dispersion (1)
params.hndsp = 1;       %     HRF negative response dispersion (1)
params.nvol = nTRs;     %     number of volumes (data points)
params.prtr = TR;       %     TR (in ms, 2000)
params.rnorm = 1;       %     Normalization (default: normalize HRF plateau to 1)
params.rcond = 0;       %     Conditions to remove (rest, 1)

data_orig=double(data);
voxGood = ~isnan(nanmean(data_orig,1));
% vthresh = 0.5;      % only let pass if GLM predicts >50% of BOLD signal & r2thresh = 1

%% get trial data
prt = BVQXfile(pname);          % PRT file labels timing of each trial
curD = data;

% load stimulus position data
mname = pname;                  % .mat file with stimulus data
mname(end-3:end) = '.mat';
load(mname,'stim');
stimPos = stim.locsRealDeg;     % stimulus position in DVA

% convolve PRT data with HRF
rtc = prt.CreateSDM(params);

%% fMRI GLM - get one beta weight for each trial
X = [rtc.RTCMatrix, ones(nTRs,1)]; % ones = constant param in RTC
% compute the beta vals for the tst set...
b = X\curD;  %inv(X'*X)*(X'*curD);

bEst = b(1:end-1,:);   % remove constant term
bEst = bEst(:,voxGood);

%% calculate how well the GLM fits the data
% (code poached from compBetas.m, JS 10/2006)
[~,x] = size(curD);
[~,tmpc]=size(X);
dm = mean(curD);            % mean of each voxel (ym)
rdm = repmat(dm,[nTRs,1]);	% repmat the mean to size of d
ssy = sum((curD-rdm).^2);   % d-rdm squared, then summed in each voxel
dfmse = nTRs-tmpc;          % degrees of freedom for determining MSE (timepoints - #predictors)
% some of these values do not need to be stored, but indexing just
% because it might be easier in the future to figure out what's going on
% (e.g. mse,ssy are not currently returned, but here we're storing these values for
% each voxel for future implementation).
for v = 1:x
   t = repmat(b(:,v)',[nTRs,1]);
   yh(:,v) = sum((t.*X)')';
   mse(v) = sum((curD(:,v) - yh(:,v)).^2)/dfmse;   %sse/dfmse
   % SE of the estimates is largely a function of design
   % matrix X, which is why we can compute the relative effeciency of 
   % a given stimulus sequence before we run the experiment...
   tmpse(:,v) = sqrt(diag(mse(v).*inv(X'*X)));             
   ssyh(v) = sum((yh(:,v)-curD(:,v)).^2);
end
r2 = 1-(ssyh./ssy);
clear tmpse ssyh yh mse dm rdm ssy dfmse

mean_r2 = mean(r2,1);

end