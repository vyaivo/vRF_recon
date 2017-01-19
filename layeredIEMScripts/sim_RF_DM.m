% sim_RF_responses.m
% given a stimulus location, radius, screen
% size, stagger offset, etc will return the response of nRFs^2 identical
% gridded uniformly-spaced bivariate cos^7 RFs (channels), 1 row per trial

% input args: stim, rfPtsX, rfPtsY, rad, fov
% - stim: contains info about each trial
% - rfPtsX: centers of basis functions in the x direction
% - rfPtsY: centers of basis functions in the y direction
% - rad: size constant of 2d cosine function

% edited 5/28/15 VAV - changed to fit with vy's code for prf project


function resp = sim_RF_DM(stim,rfPtsX,rfPtsY,rad,fov)

% rad = fwhm / rad2fwhm(1);
pow = 7; % power cos's raised to...

%% grid
% set up grid of points at which we'll sample

gridSize = 250; % sample at 250 points along x axis
xyratio = stim.maxYGridDeg / stim.maxXGridDeg;
%TODO: put this back to fov?
gridx = linspace(-stim.usedScreenSizeDeg/2, stim.usedScreenSizeDeg/2,gridSize);
gridy = linspace((-stim.usedScreenSizeDeg*xyratio)/2, (stim.usedScreenSizeDeg*xyratio)/2, gridSize*xyratio);
[x y] = meshgrid(gridx,gridy);

%% set up high-res grid of channels
% set of grid of bivariate cos^7 channels
[rfX,rfY] = meshgrid(rfPtsX, rfPtsY);
resp = zeros(numel(stim.xLocDeg),length(rfPtsX)*length(rfPtsY));
rfX = rfX(:);
rfY = rfY(:);

%% stim
resp = nan(numel(stim.xLocDeg),numel(rfX));

parfor s = 1:numel(stim.xLocDeg)
    if isnan(stim.xLocDeg(s))
        continue;
    end
    xCenter = stim.xLocDeg(s);
    yCenter = stim.yLocDeg(s);
    stimIdx = ((x-xCenter).^2 + (y-yCenter).^2).^0.5 < stim.radDeg;
% 
%     %% 
    stimX = x(stimIdx); % these are the x & y values over the area the stimulus subtends
    stimY = y(stimIdx);
    stimresp = zeros(numel(rfX),1);
    
    for ch = 1:numel(rfX)
        myr = ((stimX-rfX(ch)).^2+(stimY-rfY(ch)).^2).^0.5;
        g = ((0.5*(1 + cos(myr*pi/rad) )).^pow) .* (myr<=rad) ;
        stimresp(ch) = sum(g(:));
    end
    resp(s,:) = stimresp;

end

resp = resp./max(resp(:)); % normalize!

return