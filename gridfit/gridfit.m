% gridfit.m

function [bf_params, err, bf_fcn] = gridfit(data,grid,grid_params,amp_range)

% data should be datapts x things to fit (e.g., voxels)
% grid should be datapts x predictors
% grid_params is the params used to make grid - just used to output
% bf_params including [a b] as the last two params, and used in case
% optional fminsearch optimization step at end is used (?)
%
% uses parfor, but assumes the parallel cluster has been initialized
% already in parent funciton

if nargin < 4
    amp_range = [-inf inf];
end

use_gpu = 0;

myones = ones(size(grid,1),1);
allcoeffs = nan(size(grid,2),2,size(data,2));
allerr = nan(size(grid,2),size(data,2));


%convert things to gpuarrays
if use_gpu==1
    data = gpuArray(data);
    grid = gpuArray(grid);
    myones = gpuArray(myones);
    allcoeffs = gpuArray(allcoeffs);%nan(size(grid,2),size(data,2),2));
    allerr = gpuArray(allerr);%nan(size(grid,2),size(data,2))); % grid x vox
end


% data should be n_datapts x nvox
% grid is n_datapts x n_predictors
tic
parfor ii = 1:size(grid,2)
    allcoeffs(ii,:,:) = [grid(:,ii) myones]\data; % returns predictors(2) x nvox
    thiscoef = reshape(allcoeffs(ii,:,:),2,size(data,2));
    allerr(ii,:) = sqrt(mean((data - [grid(:,ii) myones]*thiscoef).^2,1));
end
toc

allamp = squeeze(allcoeffs(:,1,:));
allerr(allamp < amp_range(1) | allamp > amp_range(2)) = inf;
% now find the best fit and arrange those values for returning
[err,fidx] = min(allerr,[],1); % 1 x vox

if use_gpu==1
    allcoeffs = gather(allcoeffs);
    grid = gather(grid);
    err = gather(err);
end

bf_params = nan(size(data,2),size(grid_params,2)+2);
bf_fcn = nan(size(data,2),size(grid,1));    % vox x datapts
for ii = 1:length(fidx)                     % for each voxel
    bf_params(ii,:) = [grid_params(fidx(ii),:) squeeze(allcoeffs(fidx(ii),:,ii))];
    bf_fcn(ii,:)  = grid(:,fidx(ii))';
end


return