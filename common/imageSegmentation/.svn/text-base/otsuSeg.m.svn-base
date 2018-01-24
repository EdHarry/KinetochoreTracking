function mask = otsuSeg(thisimage, closureRadius)
% otsuSeg segments cell outline in fluorescence micrographs based on otsu Segmentation
% 
% Input:    thisimage                       target image for segmentation
%           (optional) closureRadius        performed image closing on the
%                                           segmentation mask with a disk
%                                           with radius = closureRadius
%
% Output:   mask                            segmented mask
%
% Last updated: May 06, 2009 by Shann-Ching Chen, LCCB
% See also: graythresh, im2bw, medianSeg, phasecontrastSeg, CustomizedSeg

if( nargin < 2)
	closureRadius = 0;
end

thisimage = nrm(thisimage,16);
level = graythresh(thisimage);
mask = im2bw(thisimage,level);

% find the largest mask
L = bwlabel(mask);
s  = regionprops(L, 'Area');
Allarea = [s.Area];
[tmp1, tmp2] = max(Allarea);
mask = (L == tmp2);

% perform closing operation
closureBrush = strel('disk',closureRadius);
mask = imclose(cast(mask,'double'),closureBrush);
mask = cast(mask,'logical');
mask = double(imfill(double(mask),'holes'));

L2 = bwlabel(mask==0);
s2  = regionprops(L2, 'Area','PixelIdxList');

% find the background with the highest average intensity
% assign the rest as forground
numberOfBlobs = size(s2, 1);
for k=1:numberOfBlobs
    PixelIdxList = s2(k).PixelIdxList;
    meanGL(k) = mean(thisimage(PixelIdxList));
end
[tmp1, tmp2] = min(meanGL);
maskbg = (L2 == tmp2);
mask = (maskbg == 0);
