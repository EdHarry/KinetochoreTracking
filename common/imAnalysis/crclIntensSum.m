function [sumTotal,avgPixInt] = crclIntensSum(im,centersYX,radii)
%CRCLINTENSSUM adds image pixel values falling within one or more circles
%
% SYNOPSIS: [sumTotal] = crclIntensSum(im,centersYX,radii)
%
% INPUT: 
%        im         : image matrix
%        centersYX  : nCircles x 2 matrix containing (row, col) indices for
%                     the center of each circle to be integrated
%        radii      : nCircles-vector containing circle radii corresponding
%                     to each circle 
%
% OUTPUT: 
%        sumTotal   : nCircles-vector containing sum of intensity for each
%                     circle. a value of NaN means the circle fell over the
%                     image boundary.
%
% This function calls crclMask, which creates a mask the size of the image
% with values between 0 (completely outside circle) and 1 (completely
% inside.  Values in between are based on how much of a given pixel on the
% edge is inside the circle.
%                     
% MATLAB VERSION (originally written on): 7.0.1.24704 (R14) Service Pack 1 Windows_NT
%
% USERNAME: kathomps
% DATE: 12-Apr-2006
%
% COMMENTS: this function works by sampling each pixel around the circle
% edge with random numbers and calculating how much of its area falls
% within the circle. pixels entirely within the circle are included (but
% not sampled).

[imL imW]=size(im); % loop thru all centers
nCircles=size(centersYX,1); % total number of centers
sumTotal=zeros(nCircles,1); % initialize vector for circle intensities to be stored

for circ=1:nCircles
    [cMask]=crclMask(imL,imW,centersYX(circ,:),radii(circ),50);
    mask=cMask.*im;
    sumTotal(circ)=sum(mask(:));
end
crclAreas=pi*radii.^2;
avgPixInt=sumTotal./crclAreas;



