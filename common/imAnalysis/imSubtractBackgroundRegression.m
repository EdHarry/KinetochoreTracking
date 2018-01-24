function [backgroundlessImage, avgBackLevel] = imSubtractBackground (varargin)
% imSubtractBackground takes an input image and subtracts the background from it. The
% image without background is returned as output. This function uses the robustfit
% algorithm to estimate a background plane and subtract it from the original image.
% This works particularly well for images with an uneven background intensity like
% phase-contrast images.
%
% SYNOPSIS       [backgroundlessImage, avgBackLevel] = imSubtractBackground (varargin)
%
% INPUT          inputImage: the original greylevel image including background
%
% OUTPUT         backgroundlessImage: the original image with the estimated background subtracted
%                avgBackLevel: the average level of the background plane
%
% DEPENDENCIES   imSubtractBackground uses { nothing }
%                                  
%                imSubtractBackground is used by { ptGetProcessedImage }
%
% Revision History
% Name                  Date            Comment
% --------------------- --------        --------------------------------------------------------
% Andre Kerstens        Apr 04          Initial release

% Do some error checking on the input
if nargin < 1
   error('The input image has to be provided. See help imSubtractBackground.');
end

inputImage = varargin{1};

% To estimate the background place we are first calculating the background plane. For this we
% need the intensity values of the 4 sides of the image and robustly fit a line through it.
% This looks as follows:
%
%                   x2 y2
% int(c1) -------------------------- int(c2)
%         |                        |
%         |                        |
%         |                        |
%      x1 |                        | x3
%      y1 |                        | y3
%         |                        |
%         |                        |
%         |                        |
% int(c4) -------------------------- int(c3)
%                   x4 y4
%

% Fetch the X input vectors for the robustfit function
x1 = (1 : size (inputImage, 1))';
x2 = (1 : size (inputImage, 2));
x3 = (1 : size (inputImage, 1))';
x4 = (1 : size (inputImage, 2));

% Fetch the Y input vectors for the robustfit function
y1 = inputImage (:,1);
y2 = inputImage (1,:);
y3 = inputImage (:,end);
y4 = inputImage (end,:);

% Calculate the coefficients of the 4 lines y=ax+b where coeff=[a,b]
coeff1 = robustfit (x1, y1);
coeff2 = robustfit (x2, y2);
coeff3 = robustfit (x3, y3);
coeff4 = robustfit (x4, y4);

% Followed by the calculation of the estimated intensity values
y1Est = coeff1(2) .* x1 + coeff1(1);
y2Est = coeff2(2) .* x2 + coeff2(1);
y3Est = coeff3(2) .* x3 + coeff3(1);
y4Est = coeff4(2) .* x4 + coeff4(1);

% These can be used to calculate average intensity values for the corners of the plane
c1 = (y1Est(1) + y2Est(1)) / 2;
c2 = (y2Est(end) + y3Est(1)) / 2;
c3 = (y3Est(end) + y4Est(end)) / 2;
c4 = (y4Est(1) + y1Est(end)) / 2;

% And create the plane itself (in matrix form)
bgPlane = [c1, c2 ; c4, c3];

% Calculate average background level
avgBackLevel = (c1 + c2 + c3 + c4) / 4;

% Ofcourse the real plane has to be subtracted from the image so we resize it
realBgPlane = imresize (bgPlane, size (inputImage), 'bilinear');

% Finally we subtract the original image with the background
backgroundlessImage = inputImage - realBgPlane;

% Let's normalize the image back to [0..1] again
imageMinimum = min (min (backgroundlessImage));
imageMaximum = max (max (backgroundlessImage));
backgroundlessImage = (backgroundlessImage - imageMinimum) / (imageMaximum - imageMinimum);

% That's it we're finished
clear x1; clear x2; clear x3; clear x4;
clear y1; clear y2; clear y3; clear y4;
clear coeff1; clear coeff2; clear coeff3; clear coeff4;
clear y1Est; clear y2Est; clear y3Est; clear y4Est;
clear c1; clear c2; clear c3; clear c4;
clear bgPlane; 
