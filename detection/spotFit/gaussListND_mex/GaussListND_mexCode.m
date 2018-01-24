function gaussList = GaussListND_mexCode(coordList,sigma,center)
%GAUSSLISTND calculates the value of a N-D Gaussian at specific pixel/voxel coordinates
%
% SYNOPSIS gaussList = GaussListND(coordList,sigma,center,intNorm,rotation)
%
% INPUT    coordList : m-by-n list of coordinates, where m is the number of
%                      coordinates and n the number of dimensions
%          sigma     : 1-by-n (or scalar): sigma of Gaussian
%          center    : (opt) 1-by-n vector of center of Gaussian.
%                      Default: zeros(1,n)
%          intNorm   : (opt) switch for how the Gaussian should be normed
%                      Default: 0
%                      0 - no norming. Max of Gaussian == 1
%                      1 - normed so that integral of infinite Gaussian = 1
%          rotation  : (opt) Equal to the number of degree you want the
%                            coordinate to be rotate for. If rotation is
%                            equal to 1, rotation will be random.
%                            Default: 0;
%                            Rotation is only supported for 2D and 3D case
%
% OUTPUT   gaussList : m-by-1 list of intensities. Intensity is the
%                      integral of the Gaussian over the pixel/voxel
%
% REMARKS  The code assumes that a pixel has the edge length 1!
%
% c: 2/05 jonas
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



[nCoords,nDims] = size(coordList);


sigma2 = repmat(sigma,[nCoords,1]);

center2 = repmat(center, [nCoords,1]);



%======================
% CALC GAUSSLIST
%======================


% instead of calculating Gauss-values for very complicated geometries, we
% make a coordinate transformation so that we can use sigma=1 in all
% dimensions



% 0.5*erfc(-(x+0.5)/sqrt(2))-0.5*erfc(-(x-0.5)/sqrt(2)) gives the integral on the
% pixel at 1 of a Gaussian with mean 0 and sigma 1


%center3 = center2(:,1);
%clear center2

% convert coordList to 0/1
coordList2 = (coordList(1:nCoords,1:nDims) - center2(1:nCoords,1:nDims))./sigma2(1:nCoords,1:nDims);
%clear coordList center3

% double coordList as preparation for erfc
%fixed bug: must divide the 0.5 by sigma - KJ
coordList2 = cat(3,coordList2-0.5./sigma2(1:nCoords,1:nDims), coordList2+0.5./sigma2(1:nCoords,1:nDims));

% calculate gaussList
%Jonas was missing the minus sign in erfc. I corrected that - KJ
gaussList = diff(0.5 * erfc(-coordList2/sqrt(2)),1,3);
gaussList = prod(gaussList,2);

% norm gaussList
gaussList = gaussList*((2*pi)^(0.5*nDims)*prod(sigma2(1,:)));

