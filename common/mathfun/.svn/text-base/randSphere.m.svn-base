function coordinates = randSphere(nPoints, nDims, radius)
%RANDSPHERE returns uniformly distributed random points within a sphere
%
% SYNOPSIS coordinates = randSphere(nPoints, nDims, radius)
%
% INPUT    nPoints : number of random points
%          nDims   : (opt) number of dimensions of the sphere (2 = circle)
%                       Default: 3
%          radius  : (opt) radius of the sphere. If length(radius)==nDims,
%                       the radius can be set individually for each
%                       dimension. Default: 1
%
% OUTPUT   coordinates : nPoints-by-nDims array with point coordinates
%
% c: jonas, 11/05
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%===================
% TEST INPUT
%===================

% set defaults
def_nDims = 3;
def_radius = 1;
% number of potential coordinates the program will produce in each
% iteration
pointMultiplicator = 3; 

% check nargin
if nargin == 0
    error('please specify at least the input argument ''nPoints''')
end

if nargin < 2 || isempty(nDims)
    nDims = def_nDims;
end

if nargin < 3 || isempty(radius)
    radius = repmat(def_radius,1,nDims);
else
    % check for correct number of elements
    radius = radius(:)';
    switch length(radius)
        case 1 % expand to all the same
            radius = repmat(radius,1,nDims);
        case nDims
            % perfect
        otherwise
            error(['radius needs to either be a '...
                'scalar or have nDims (=%i) elements"'],nDims);
    end
end

%================

%================
% CREATE POINTS
%================

% to make correctly distributed points, we generate 3xnPoints coordinates,
% and remove all that are outside the circle. If we don't have enough
% points, we just do another 3xn points, until we're done
% First, coordinates are generated for a cube with center at 0 whose sides
% have length 2. Then, the values are scaled by the radii (division by
% normed radius). Every point that is further than 1 away from the center
% is discarded from the original list. The remaining points are then
% multiplied by the radius.

done = 0;
coordinates = zeros(nPoints,nDims);
lastCoordinate = 0;
maxRadius = max(radius);
radiusRatio = repmat(...
    radius./maxRadius, ...
    nPoints * pointMultiplicator, 1);

while ~done
    % generate potential coords
    potentialCoordinates = 2* rand(nPoints * pointMultiplicator,nDims) - 1;
    
    % get length of coords that have been scaled by the inverse of the
    % radius
    distanceFromOrigin = ...
        sqrt(sum((potentialCoordinates ./ radiusRatio).^2,2));
    
    % find good coords
    goodCoordIdx = find(distanceFromOrigin <= 1);
    
    % write into output
    nGoodCoords = length(goodCoordIdx);
    % make sure we don't write out too much
    maxIdx = nPoints - lastCoordinate;
    if maxIdx < nGoodCoords
        % too many good coords: only use the first maxIdx
        nGoodCoords = maxIdx;
        goodCoordIdx = goodCoordIdx(1:maxIdx);
    end
    
    coordinates(lastCoordinate + 1:lastCoordinate + nGoodCoords, :) = ...
        potentialCoordinates(goodCoordIdx,:);
    
    % update lastCoordinate
    lastCoordinate = lastCoordinate + nGoodCoords;
    
    % quit if done
    if lastCoordinate == nPoints
        done = 1;
    end
    
end % while ~done

% adjust radius
coordinates = coordinates .* maxRadius;