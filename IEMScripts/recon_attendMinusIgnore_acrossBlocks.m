function [fdat] = recon_attendMinusIgnore_acrossBlocks(sfit,atlocs,condlocs)
% collapses center fit values (x,y) into a single distance from attention
% locus parameter. then finds the difference between all parameters
% (distance from attention locus, size, amplitude, baseline) across runs of
% the attention task (e.g., L1 when attending left - L1 when attending
% right).
% function useful for plotting & stats of reconstructions
% INPUTS
% sfit: niterations x nconditions x npositions x nparameters (fit params,
% e.g. 100 x 2 x 51 x 5) or can be 2 x 51 x 5
% atlocs: 2 (L/R) x 2 (x/y coordinates of attention loci)
% condlocs: 2 (L/R) x npairedlocs (e.g. 2 x 24). specifies which positions
% are paired across conditions (e.g., condlocs(:,1) is the first paired set
% of positions across hemifields

if numel(size(sfit)) == 4

    niters = size(sfit,1);
    nc = size(sfit,2);
    % initialize output matrix
    fdat = nan(niters,nc,size(condlocs,2),size(sfit,4)-1);

    for c = 1:nc
        ic = setxor(c,1:1:nc);
        % first set the dist from attend param
        x1 = squeeze(sfit(:,c,condlocs(c,:),1));
        y1 = squeeze(sfit(:,c,condlocs(c,:),2));
        x2 = squeeze(sfit(:,ic,condlocs(c,:),1));
        y2 = squeeze(sfit(:,ic,condlocs(c,:),2));
        atdist = sqrt( (x1-atlocs(c,1)).^2 + (y1-atlocs(c,2)).^2 );
        igdist = sqrt( (x2-atlocs(c,1)).^2 + (y2-atlocs(c,2)).^2 );
        fdat(:,c,:,1) = atdist - igdist;
        % now fill in the other params
        atpar = sfit(:,c,condlocs(c,:),3:end);
        igpar = sfit(:,ic,condlocs(c,:),3:end);
        fdat(:,c,:,2:end) = atpar - igpar;
    end
    
else
    nc = size(sfit,1);
    % initialize output matrix
    fdat = nan(nc,size(condlocs,2),size(sfit,3)-1);

    for c = 1:nc
        ic = setxor(c,1:1:nc);
        % first set the dist from attend param
        x1 = squeeze(sfit(c,condlocs(c,:),1));
        y1 = squeeze(sfit(c,condlocs(c,:),2));
        x2 = squeeze(sfit(ic,condlocs(c,:),1));
        y2 = squeeze(sfit(ic,condlocs(c,:),2));
        atdist = sqrt( (x1-atlocs(c,1)).^2 + (y1-atlocs(c,2)).^2 );
        igdist = sqrt( (x2-atlocs(c,1)).^2 + (y2-atlocs(c,2)).^2 );
        fdat(c,:,1) = atdist - igdist;
        % now fill in the other params
        atpar = sfit(c,condlocs(c,:),3:end);
        igpar = sfit(ic,condlocs(c,:),3:end);
        fdat(c,:,2:end) = atpar - igpar;
    end
    
end