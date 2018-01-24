function [corfac] = circumferenceCorrFactor_MaskEU(xx,yy,rr,msx,msy,mask)
% circumference correction calculates a vector containing the correction factor
% (for edge correction in Ripley's k-function) for values of rr
% circumference correction: fraction of circumference of circle centered at
% point P=(xx,yy) with radius rr (inside rectangular image) falling into the 
% area of interest defined by the mask image 
%
% SYNOPSIS
% [corfac]=circumferenceCorrFactor_MaskEU(xx,yy,rr,msx,msy,mask)
%
% INPUT:    xx  = x-position
%           yy  = y-position
%           rr  = distance vector
%           msx = image size in x
%           msy = image size in y
%           mask=region-of-interest image mask
%           
%           NOTE is mask is a single value rather than a matrix, it's
%           assumed that the distance transform is performed outside the
%           function and that mask alerady represents the euclidian
%           distance value for this point
%
% OUTPUT:   corfac = correction factor values (vector of same length as rr)
%
% Notes: - this function assumes that rr is a vector
%        - this version of the function allows rr to be larger than ms/2
%       
% Dinah Loerke, Feb 05, 2008
%
% NOTE: to increase speed, this version of the function uses the following
% approximations: the area of interest is considered to be a roundish blob;
% locally for each position inside the blob, the circumference ratio is
% considered to be that of a circle cut off by a straight edge at distance
% ed, which is the euclidian distance from the mask egde. This approach has
% the advantage that the euclidian distance can be calculated quickly with
% bwmask, and that the straight edge cutoff of a circle can be calculated
% analytically from geometrical considerations


% due to the grid-based/sample approach to the correction factor problem,
% it can happen occasionally that the chosen coordinates xx,yy lie OUTSIDE
% the area mask; there are two possibilities for this case, on the one
% hand, cofac values could be set to nan, on the other hand, cofac values
% could be set to the values for the CLOSEST point to xx,yy inside the area

px = round(xx);
py = round(yy);

% if mask is single value, use this for euclidian distance, else calculate
% the distance matrix here
if length(mask(:))==1
    ed = mask;
else
    % the dimensions of the mask have to be the same as msx, msy
    [mx,my] = size(mask);
    if (mx~=msx) | (my~=msy)
        error('image and mask size don''t agree');
    end
    maskED = bwdist(mask);
    ed = maskED(px,py);
end


%correction factor is initialized
corfac=1+0*rr;
rl = length(rr);

% minimum position for which circle radius is larger than euclidian
% distance; if the circle radius is smaller, nothing is cut off and the
% correction factor equals one
mpos = find(rr>ed);

if ~isempty(mpos)
    
    for ri = min(mpos):rl
        % as in geometrical factor calculation
        % current distance
        currDist = rr(ri);
        circRatio = 0.5 + 0.5 * ( asin(ed/currDist) / (pi/2) );

        corfac(ri) = circRatio;

    end % of for ri
    
end % of if

    

end  % of function


