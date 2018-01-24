function mask = phasecontrastSeg(im, closureRadius)
% phasecontrastSeg segments cell outline in phase contrast micrographs
% Use |edge| and the Sobel|log|canny operator to calculate the threshold
% value. Tune the threshold value and use |edge| again to obtain a
% binary mask that contains the segmented cell. This method is modifying
% ipexbatchDetectCells.m in the image processing toolbox
% 
% Input:    im                              target image for segmentation
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

im1 = edge(im,'log');
im2 = edge(im,'sobel');
im3 = edge(im,'canny');
%im4 = tz_imstdedge(im)>0;
%im5 = ipexSeg(im);

%BW = im1 | im2 | im3 | im4 | im5;
BW = im1 | im2 | im3;

se90 = strel('line', 3, 90);
se0 = strel('line', 3, 0);
BWdilate = imdilate(BW, [se90 se0]);
BWdilate([1,end],:) = 0;  BWdilate(:,[1,end]) = 0;
BWnobord = imclearborder(BWdilate, 4);
BWnobord(1,:) = BWnobord(2,:);  
BWnobord(:,1) = BWnobord(:,2);
BWnobord(end,:) = BWnobord(end-1,:);  
BWnobord(:,end) = BWnobord(:,end-1);

BWopen = bwareaopen(BWnobord,200);
BWclose = bwmorph(BWopen,'close');
BWclose(1,:) = BWclose(2,:);  
BWclose(:,1) = BWclose(:,2);
BWclose(end,:) = BWclose(end-1,:);  
BWclose(:,end) = BWclose(:,end-1);

L = bwlabel(BWclose);
s  = regionprops(L, 'Area');
Allarea = [s.Area];
[tmp1, tmp2] = max(Allarea);
mask = (L == tmp2);

idx = find(mask(1,:)==1);    mask(1,min(idx):max(idx)) = 1;
idx = find(mask(:,1)==1);    mask(min(idx):max(idx),1) = 1;
idx = find(mask(end,:)==1);  mask(end,min(idx):max(idx)) = 1;
idx = find(mask(:,end)==1);  mask(min(idx):max(idx),end) = 1;

closureBrush = strel('disk',closureRadius);
mask = imclose(cast(mask,'double'),closureBrush);
mask = cast(mask,'logical');

mask = bwfill(mask, 'holes', 4);

L2 = bwlabel(mask==0);
s2  = regionprops(L2, 'Area','PixelIdxList');

numberOfBlobs = size(s2, 1);
for k=1:numberOfBlobs
    PixelIdxList = s2(k).PixelIdxList;
    meanGL(k) = mean(mask(PixelIdxList));
end
[tmp1, tmp2] = min(meanGL);
maskbg = (L2 == tmp2);
mask = (maskbg == 0);
