function basis_set = build_basis_pts(centersX,centersY,basis_size,xx,yy,varargin)

% Builds basis set given just a list of basis fcn centers, basis fcn sizes, 
% and an x,y grid along which to evaluate the 2dcos (defined by make2dcos).

% This is a wrapper around /gridfit/make_grid.m, which define surfaces
% specified by an input function (here, make2dcos) in a specific grid of
% points (here, meshgrid(centersX,centersY)).

% Tommy Sprague (TCS), 10/24/14
% Altered VAV 1/19/15 to accomodate a RECTANGULAR basis set

if nargin > 5
    parflag = varargin{1};
elseif nargin == 5
    % unless you turn it off, this will make each basis fcn on a
    % separate parallel processor. good when trying to make a huge number
    % of functions for a grid search to fit data.
    parflag = 1;
end

if length(basis_size)==1
    basis_size = basis_size*ones(length(centersX),1);
end

basis_set = make_grid(@make2dcos_grid,[xx yy],[centersX centersY basis_size],parflag);


return