function mask = medianSeg(im, closureRadius)
% medianSeg segments cell outline in fluorescence micrographs based on
% median assignment followed by otsu Segmentation 
% 
% Input:    im                              target image for segmentation
%           (optional) closureRadius        performed image closing on the
%                                           segmentation mask with a disk
%                                           with radius = closureRadius
%
% Output:   mask                            segmented mask
%
% Last updated: May 06, 2009 by Shann-Ching Chen, LCCB
% See also: otsuSeg, phasecontrastSeg, CustomizedSeg

if( nargin < 2)
	closureRadius = 0;
end

% apply guassian smoothing
img_org = im;
h = fspecial('gaussian',[7 7], 2);
img_org = imfilter(img_org,h,'replicate');

% supressed the image intensity : for those pixels with intensities larger
% than the median, assign the intensity as median
Q(2) = median(im(:));
img_org(find(im>Q(2))) = Q(2);

thismask = nrm(img_org,16);
level = graythresh(thismask);

if level > 0
    mask = im2bw(thismask,level);
    L = bwlabel(mask);
    s  = regionprops(L, 'Area');
    Allarea = [s.Area];
    [tmp1, tmp2] = max(Allarea);
    mask = (L == tmp2);
    
    closureBrush = strel('disk',closureRadius);
    mask = imclose(cast(mask,'double'),closureBrush);
    mask = cast(mask,'logical');

    mask = double(imfill(double(mask),'holes'));
    
    if length(unique(mask)) == 1
        mask = thismask > median(double(thismask(:)));
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
    end
else
    % perform closing operation    
    mask = thismask > median(thismask(:));
    L = bwlabel(mask);
    s  = regionprops(L, 'Area');
    Allarea = [s.Area];
    [tmp1, tmp2] = max(Allarea);
    mask = (L == tmp2);

    closureBrush = strel('disk',closureRadius);
    mask = imclose(cast(mask,'double'),closureBrush);
    mask = cast(mask,'logical');

    mask = double(imfill(double(mask),'holes'));
end

L2 = bwlabel(mask==0);
s2  = regionprops(L2, 'Area','PixelIdxList');

% find the background with the highest average intensity
% assign the rest as forground
numberOfBlobs = size(s2, 1);
for k=1:numberOfBlobs
    PixelIdxList = s2(k).PixelIdxList;
    meanGL(k) = mean(thismask(PixelIdxList));
end
[tmp1, tmp2] = min(meanGL);
maskbg = (L2 == tmp2);
mask = (maskbg == 0);

