function [minThresh, J] = imMinimumThreshold (image, maxlevel)
% imMinimumThreshold computes the minimum threshold value needed to effectively
% segment that image.  The assumption is made that the objects and pixel grey
% values are normally distributed. This algorithm was first designed by J. Kittler
% and J. Illingworth.
%
% SYNOPSIS [J] = imMinimumThreshold (image, maxlevel)
%
% INPUT  image:     The input image in normalized form (0..1)
%        maxlevel:  True bit depth of the image (eg 255 for 8 bit images)
%
% OUTPUT J:         the values of the criterium function in vector form
%        minThresh: the minimum threshold value found
%
% Revision History
% Name                  Date            Comment
% --------------------- --------        --------------------------------------------------------
% Andre Kerstens        May 04          Initial release


if (nargin < 2)
   error ('This function needs two parameters: a normalized image and the maximum greylevel value (e.g. 4095 for a 12 bit image). See help imMinimumThreshold.');
end;

% De-normalize the image to grey level values
imageGrey = image.*maxlevel;

% Make sure the threshold value actually makes sense for the image. It
% should be above the min greylevel and under max greylevel of the image
threshMin = min (min (imageGrey)) + 3;
threshMax = max (max (imageGrey)) - 3;

% Calculate the histogram of the input image
[h, loc] = imhist (image, maxlevel+1);

% generate a vector of the greyvalues in the image
grey = 0:1:maxlevel;

% Initialize the index for J
indJ = 1;

% Calculate the values for the criterion function J
%for iCount = threshMin : 1 : threshMax
for iCount = 1 : 1 : maxlevel+1
    
   histBelowThresh = h (1:iCount);
   histAboveThresh = h (iCount+1:maxlevel+1);

   % Take the parts of the greylevels that are below and under the threshold
   imageBelowThresh = grey (1:iCount)';
   imageAboveThresh = grey (iCount+1:maxlevel+1)';

   % Calculate P1(T) and P2(T)
   p1 = sum (histBelowThresh);
   p2 = sum (histAboveThresh);

   % Calculate mu1(T), mu2(T), sigma1(T) and sigma2(T)
   if p1 ~= 0
      mu1 = sum (histBelowThresh .* imageBelowThresh) / p1;
      sigma1 = sum ((imageBelowThresh - mu1).^2 .* histBelowThresh) / p1;
   else
      mu1 = 0;
      sigma1 = 0;
   end
   if p2 ~= 0
      mu2 = sum (histAboveThresh .* imageAboveThresh) / p2;
      sigma2 = sum ((imageAboveThresh - mu2).^2 .* histAboveThresh) / p2;
   else
      mu2 = 0;
      sigma2 = 0;
   end

   % Calculate the criterion function J
   if (p1 ~= 0) & (p2 ~= 0) & (sigma1 ~= 0) & (sigma2 ~= 0)
      J(indJ) = 1 + 2 * (p1 * log(sigma1) + p2 * log(sigma2)) - 2 * (p1 * log(p1) + p2 * log(p2));
   else
      J(indJ) = 0;
   end
   
   indJ = indJ + 1;
end

% Get the minimum value of J
[minJ, minIndex] = min (J);
minThresh = grey (minIndex) / maxlevel;
