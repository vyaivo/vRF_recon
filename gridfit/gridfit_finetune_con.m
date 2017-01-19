function [bf_ft,bf_err,bf_fcn,ex_flag] = gridfit_finetune_con(data,fitfcn,...
    bf_grid,evalpts,gridstep,constr,parflag,shutup)

% refines fits computed with gridfit.m using fminsearch and fitfcn, which
% is used here to generate an error function (root mean squared error)
%
% assumes bf_grid is output from gridfit.m - last 2 columns are amplitude
% and baseline, respectively
%
% data formatted like gridfit.m:
% - data: n_datapts x n_things_to_fit (voxels)
% - fitfcn: function handle to a normalized-amplitude/baseline function
%   prototype used for GLM fitting in gridfit fitting (base = 0, amp = 1)
% - bf_grid: the best-fit from gridfit.m (its return argument bf_params)
% - evalpts: same as gridfit.m, the points at which fitfcn is evaluated for
%   a given set of parameters (e.g, x/y coords) - n_datapts x n_coords
% - gridstep: a vector (1 x n_fit_params - 2) specifying how much the fit
%   can deviate from the best fit in bf_grid. does not include amplitude &
%   baseline, which have separate constraints (passed in as constr)
% - constr: a 2 x 2 matrix of [lower bounds; upper bounds] on
%   the amplitude & baseline of the fits.
% - parflag: 1 if you want parallel processing, 0 if not
% - shutup: 0 if you want output, 1 if not

% TCS 11/2014; VAV major edits 2/2015

if nargin < 8    
    shutup = 0;
end
if nargin < 7
    parflag = 1;
end
if nargin < 6
    constr = [0,-5; 5,5];
end

% Make the LB & UB list (it's just bestfit +/- gridstep)
LB = bsxfun(@minus,bf_grid(:,1:length(gridstep)),gridstep);
UB = bsxfun(@plus,bf_grid(:,1:length(gridstep)),gridstep);

% initialize all variables
bf_ft = nan(size(bf_grid));
bf_err = nan(size(bf_grid,1),1);
bf_fcn = nan(size(data));
ex_flag = nan(size(bf_grid,1),1);

amp_idx = size(bf_grid,2)-1;
base_idx = size(bf_grid,2);


% make fmincon shut up
options = optimset('Display','off','Algorithm','active-set');
options.Display = 'off';


if parflag
    parfor vv = 1:size(data,2)

        if ~shutup && mod(vv,10)==0
            fprintf('Fit %i\n',vv);
        end

        d = data(:,vv);

        % here's the error function (RMSE)
        % p(end-1) is amplitude & p(end) is baseline. so here it's just
        % sqrt(mean( ((amp*fitfunction + baseline) - data)).^2 ))
        err_fcn = @(p) sqrt(mean( (   (p(end-1)*fitfcn(evalpts,p)+p(end)) - d).^2 ) );

        [bf_ft(vv,:),bf_err(vv),ex_flag(vv)] = fmincon(err_fcn,bf_grid(vv,:),...
            [],[],[],[],[LB(vv,:) constr(1,:)],[UB(vv,:) constr(2,:)],[],options);
        this_fits = bf_ft(vv,:);
        this_amp = this_fits(amp_idx);
        this_base = this_fits(base_idx);
        bf_fcn(:,vv) = this_amp*fitfcn(evalpts,bf_ft(vv,:)) + this_base;

    end
else
    for vv = 1:size(data,2)

        if ~shutup && mod(vv,10)==0
            fprintf('Fit %i\n',vv);
        end

        d = data(:,vv);

        % here's the error function (RMSE)
        err_fcn = @(p) sqrt(mean( (   (p(end-1)*fitfcn(evalpts,p)+p(end)) - d).^2 ) );

        [bf_ft(vv,:),bf_err(vv),ex_flag(vv)] = fmincon(err_fcn,bf_grid(vv,:),...
            [],[],[],[],[LB(vv,:) constr(1,:)],[UB(vv,:) constr(2,:)],[],options);
        this_fits = bf_ft(vv,:);
        this_amp = this_fits(amp_idx);
        this_base = this_fits(base_idx);
        bf_fcn(:,vv) = this_amp*fitfcn(evalpts,bf_ft(vv,:)) + this_base;

    end
end


return