function [gridCoordsYX]=gridPtCoordsYX(roiL, roiW, intY, intX)
% DESCRIPTION: returns YX-coordinates for a grid centered in a box of given
%              dimensions
%
% SYNOPSIS: [gridCoordsYX]=gridPtCoordsYX(roiL, roiW, intY, intX)
%
% INPUT:    roiL/W: length and width of region on which to create grid 
%           intY/X: grid spacing in row (Y) and column (X) directions
%
% OUTPUT:   gridCoordsYX: n x 2 matrix containing YX-coordinates of the
%                         grid points
%
% MATLAB VERSION (originally written on): 7.2.0.232 (R2006a) Windows_NT
% USERNAME: kathomps DATE: 31-Oct-2006
%

% calculate how many grid spacings fit in row and column directions
% (will end up with n spaces between n+1 grid points)
nGridSpacingsR = floor((roiL-1)/intY);
nGridSpacingsC = floor((roiW-1)/intX);

% find the first point so that the grid is centered
firstR = ceil((roiL-nGridSpacingsR*intY)/2);
firstC = ceil((roiW-nGridSpacingsC*intX)/2);

% last point in row and column directions
lastR = firstR + nGridSpacingsR*intY;
lastC = firstC + nGridSpacingsC*intX;

% make nx2 matrix gridCoordsYX of grid coordinates, where n = total num of
% grid positions.  gridCoordsYX is in (y,x) form.
ptsR = linspace(firstR,lastR,nGridSpacingsR+1);
ptsC = linspace(firstC,lastC,nGridSpacingsC+1);
[C,R] = meshgrid(ptsC,ptsR);
Y = R(:);
X = C(:);
gridCoordsYX = [Y X];