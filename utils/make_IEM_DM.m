function [X,basis_set,fov,res] = make_IEM_DM(stim,chanX,chanY,fwhm,fov,res)
% X = make_IEM_DM(stim,chanX,chanY,fwhm,fov,res)
% 
% This is a function to generate the IEM design matrix based on the known
% parameters for the basis functions (which are assumed to be 2d cosines,
% as defined in make2dcos.m) and the trial-by-trial stimulus locations.
% These are taken from a .mat output file from each run of the experimental
% paradigm.

% Inputs:
% stim is a structure from the behavioral file
% chanX & chanY are the x,y centers of the basis functions
% FOV will be defined by the basis set

if nargin < 5
    fov(1) = 2*(max(chanX) + fwhm);
    fov(2) = 2*(max(chanY) + fwhm);
end

% if no resolution defined, make sure square-ish pixels
if nargin < 6
    if length(chanX) < length(chanY)
        res = [101 round((length(chanY)/length(chanX)) * 101)];
        if mod(res(2),2) == 0
            res(2) = res(2)+1;
        end
    elseif length(chanY) < length(chanX)
        res = [round((length(chanX)/length(chanY)) * 101) 101];
        if mod(res(1),2) == 0
            res(1) = res(1)+1;
        end
    else
        res = [101 101];
    end
end

filt_size =  (1/rad2fwhm(1)) * fwhm;    % convert fwhm to size constant

% pixel grid over which stimulus mask & basis set is defined
gridxPts = linspace(-fov(1)/2,fov(1)/2,res(1));
gridyPts = linspace(-fov(2)/2,fov(2)/2,res(2));
[gridx, gridy] = meshgrid(gridxPts,gridyPts);

[rfx, rfy] = meshgrid(chanX,chanY);
% build our basis set (nChannels x nPixels)
basis_set = build_basis_pts(rfx(:),rfy(:),filt_size,gridx(:),gridy(:),0);

% now, we need to make a stimulus filter for each trial
% stim.xLocDeg, yLocDeg contains all center stim position in screen coords
% stim.radDeg is size of stimulus
stim_mask = nan(size(stim.locsRealDeg,1),res(1)*res(2));

for ss = 1:size(stim.locsRealDeg,1)
    % make a circular mask for each checkerboard stimulus
    thisr = sqrt((gridx-stim.locsRealDeg(ss,1)).^2 + (gridy-stim.locsRealDeg(ss,2)).^2);
    thisstim = thisr <= stim.radDeg;
    % check by imagesc(thisstim)
    stim_mask(ss,:) = reshape(thisstim,1,numel(thisstim));
    clear thisstim thisr;
end

% stim_mask is nTrials x nPixels
% we want nTrials x nChannels

X = stim_mask * basis_set;

% intentionally not normalizing X here, it should be normalized across all
% runs, and this will work on a single run...

return