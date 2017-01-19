 function hh = quiverc_away_toward_VAV(varargin)
 
% Modified version of Quiver to plots velocity vectors as arrows 
% with components (u,v) at the points (x,y).
% VAV: Further modified to plot arrows pointing toward some center point as
% colormap 1 (default 'winter') and away from some center point as colormap 2
% (default 'autumn').

% Vy Vo 10/3/2015
% Bertrand Dano 3-3-03
% Copyright 1984-2002 The MathWorks, Inc. 

%QUIVERC Quiver color plot.
%   QUIVERC(X,Y,U,V) plots velocity vectors as arrows with components (u,v)
%   at the points (x,y).  The matrices X,Y,U,V must all be the same size
%   and contain corresponding position and velocity components (X and Y
%   can also be vectors to specify a uniform grid).  QUIVER automatically
%   scales the arrows to fit within the grid.
%
%   QUIVERC(U,V) plots velocity vectors at equally spaced points in
%   the x-y plane.
%
%   QUIVERC(U,V,C) or QUIVER(X,Y,U,V,C) figures out whether the vector
%   direction toward or away from a center point C (default [0,0]).
%
%   QUIVERC(...,LINESPEC) uses the plot linestyle specified for
%   the velocity vectors.  Any marker in LINESPEC is drawn at the base
%   instead of an arrow on the tip.  Use a marker of '.' to specify
%   no marker at all.  See PLOT for other possibilities.
%
%   QUIVERC(...,'filled') fills any markers specified.
%
%   H = QUIVERC(...) returns a vector of line handles.
%
%   Example:
%      [x,y] = meshgrid(-2:.2:2,-1:.15:1);
%      z = x .* exp(-x.^2 - y.^2); [px,py] = gradient(z,.2,.15);
%      contour(x,y,z), hold on
%      quiverc(x,y,px,py), hold off, axis image
%
%   See also FEATHER, QUIVER3, PLOT. 
%   Clay M. Thompson 3-3-94
%   Copyright 1984-2002 The MathWorks, Inc. 
%   $Revision: 5.21 $  $Date: 2002/06/05 20:05:16 $ 
%-------------------------------------------------------------

% set(gca, 'color', 'blue');
% Arrow head parameters
alpha = 0.4; % Size of arrow head relative to the length of the vector
% alpha = 0.4;
beta = 0.6;  % Width of the base of the arrow head relative to the length
autoscale = 0; % Autoscale if ~= 0 then scale by this.
ctr = [0 0];
plotarrows = 1; % Plot arrows
sym = '';

filled = 0;
ls = '-';
ms = '';
col = '';
lw=2;           % line width

nin = nargin;
% Parse the string inputs
while isstr(varargin{nin}),
  vv = varargin{nin};
  if ~isempty(vv) & strcmp(lower(vv(1)),'f')
    filled = 1;
    nin = nin-1;
  else
    [l,c,m,msg] = colstyle(vv);
    if ~isempty(msg), 
      error(sprintf('Unknown option "%s".',vv));
    end
    if ~isempty(l), ls = l; end
    if ~isempty(c), col = c; end
    if ~isempty(m), ms = m; plotarrows = 0; end
    if isequal(m,'.'), ms = ''; end % Don't plot '.'
    nin = nin-1;
  end
end


% error(nargchk(2,5,nin));

% Check numeric input arguments
if nin<4, % quiver(u,v) or quiver(u,v,s)
  [msg,x,y,u,v] = xyzchk(varargin{1:2});
else
  [msg,x,y,u,v] = xyzchk(varargin{1:4});
end
if ~isempty(msg), error(msg); end

if nin==3 | nin==5, % quiver(u,v,s) or quiver(x,y,u,v,s)
%   autoscale = varargin{nin};
    ctr = varargin{nin};
end

% Scalar expand u,v
if prod(size(u))==1, u = u(ones(size(x))); end
if prod(size(v))==1, v = v(ones(size(u))); end

if autoscale,
  % Base autoscale value on average spacing in the x and y
  % directions.  Estimate number of points in each direction as
  % either the size of the input arrays or the effective square
  % spacing if x and y are vectors.
  if min(size(x))==1, n=sqrt(prod(size(x))); m=n; else [m,n]=size(x); end
  delx = diff([min(x(:)) max(x(:))])/n;
  dely = diff([min(y(:)) max(y(:))])/m;
  len = sqrt((u.^2 + v.^2)/(delx.^2 + dely.^2));
  autoscale = autoscale*0.9 / max(len(:));
  u = u*autoscale; v = v*autoscale;
end

%----------------------------------------------
% VAV modified in this chunk, 10/3/2015

% Define colormap
vr=sqrt(u.^2+v.^2);
dist1 = sqrt( (x-ctr(1)).^2 + (y-ctr(2)).^2 );
dist2 = sqrt( ((x+u)-ctr(1)).^2 + ((y+v)-ctr(2)).^2 );
% toflag: 1 if moving toward, -1 if moving away
toflag = sign(dist1-dist2);

% I want to color inward pointing vectors separately from outward pointing
% vectors.
cmap1 = colormap('winter');
cmap2 = colormap('autumn');
CC = cat(3,cmap1,cmap2);
% CC = repmat([0 0 0],64,1,2);

vrn=round(vr/max(vr)*64);
% CC=colormap;
ax = newplot;
next = lower(get(ax,'NextPlot'));
hold_state = ishold;

%----------------------------------------------
% Make velocity vectors and plot them

x = x(:).';y = y(:).';
u = u(:).';v = v(:).';
vrn=vrn(:).';
uu = [x; x+u; repmat(NaN,size(u))];
vv = [y; y+v; repmat(NaN,size(u))];
vrn1= [vrn; repmat(NaN,size(u)); repmat(NaN,size(u))];
toflag1 = [toflag; toflag; toflag];
% vrn1 = repmat(vrn,3,1);
% vrn1 = vrn1(:);
% toflag1 = repmat(toflag',3,1);
% toflag1 = toflag1(:);

uui=uu(:);  vvi=vv(:);  vrn1=vrn1(:); imax=size(uui); toflag1=toflag1(:);
hold on

for i = 1:3:imax-1
    ii=int8(round(vrn1(i)));
    % default colormap
    ci = 2;
    if ii==0; ii=1; end
%     c1= CC(ii,1);    c2= CC(ii,2);    c3= CC(ii,3);
    if toflag1(i)==1; ci=1; end
    h1(i) = plot(uui(i:i+1),vvi(i:i+1),'LineWidth',lw,'Color',CC(ii,:,ci));
end

%----------------------------------------------
% Make arrow heads and plot them
if plotarrows,
 
  hu = [x+u-alpha*(u+beta*(v+eps)); x+u; ...
        x+u-alpha*(u-beta*(v+eps)); repmat(NaN,size(u))];
  hv = [y+v-alpha*(v-beta*(u+eps)); y+v; ...
        y+v-alpha*(v+beta*(u+eps)); repmat(NaN,size(v))];
  vrn2= [vrn;vrn;vrn;vrn];
  toflag2 = [toflag;toflag;toflag;toflag];
%   vrn2 = repmat(vrn,4,1);
%   vrn2 = vrn2(:);
%   toflag2 = repmat(toflag',4,1);
%   toflag2 = toflag2(:);

  hui=hu(:);  hvi=hv(:);  vrn2=vrn2(:); imax=size(hui);

 for i = 1:imax-1
    ii=int8(round(vrn2(i)));
    if ii==0; ii=1; end
%     c1= CC(ii,1);    c2= CC(ii,2);    c3= CC(ii,3);
    ci = 2; if toflag2(i)==1; ci=1; end
    h2(i) = plot(hui(i:i+1),hvi(i:i+1),'LineWidth',lw,'Color',CC(ii,:,ci));
 end

else
  h2 = [];
end
%----------------------------------------------

if ~isempty(ms), % Plot marker on base
  hu = x; hv = y;
  hold on
  h3 = plot(hu(:),hv(:),[col ms]);
  if filled, set(h3,'markerfacecolor',get(h1,'color')); end
else
  h3 = [];
end

if ~hold_state, hold off, view(2); set(ax,'NextPlot',next); end

% if nargout>0, hh = [h1;h2;h3]; end
hh = struct('lines',[],'arrowheads',[],'markers',[]);
hh.lines = h1;
hh.arrowheads = h2;
hh.markers = h3;
% set(gca, 'color', [0 0 0],'Xcolor','w','Ycolor','w');
% set(gcf, 'color', [0 0 0]);
%set(gcf, 'InvertHardCopy', 'off');
