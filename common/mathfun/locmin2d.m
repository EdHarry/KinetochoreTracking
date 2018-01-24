function fImg=locmin2d(img,mask)
%LOCMIN2D searches for local minima in an image
%
%    SYNOPSIS fImg=(img,mask)
%
%    INPUT    img    image matrix
%             mask   [m n] defines the operator window dimensions
%
%    OUTPUT   fImg   map with all the local minima (zeros elsewhere);
%                    the non-zero values contain the original value 
%                    of the image at that place
%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% PARAMETER CHECK
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin<2
   error('Please define all parameters');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% DEFINITIONS
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% make sure the mask elements are odd numbers (only then, the local max operator is 
% properly defined)
indx = find(~mod(mask,2));
mask(indx) = mask(indx) + 1;

% apply a min filter
fImg = ordfilt2(img,1,ones(mask));
fImg2 = ordfilt2(img,2,ones(mask));
fImg(find(fImg2==fImg))=0;

% take only those positions where the min filter and the original image value
% are equal -> this is a local minimum
fImg(~(fImg == img)) = 0;

% set image border to zero
auxM = zeros(size(img));
auxM(fix(mask(1)/2)+1:end-fix(mask(1)/2),fix(mask(2)/2)+1:end-fix(mask(2)/2)) = ...
    fImg(fix(mask(1)/2)+1:end-fix(mask(1)/2),fix(mask(2)/2)+1:end-fix(mask(2)/2));
fImg=auxM;



