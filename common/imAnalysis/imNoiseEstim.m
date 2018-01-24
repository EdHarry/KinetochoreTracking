function [nse,thresh,grad]=imNoiseEstim(img,c)
%IMNOISEESTIM returns an estimate of the image noise
% The algorithms is described in Voorhees and Poggio, 
%                                Darpa Image Understanding Workshop 892-899, 1987.
%
% SYNOPSIS [nse,thresh,grad]=imNoiseEstim(img,c)
%
% INPUT img : greyvalue image (uint8 or double accepted)
%       c   : confidence probability on which the threshold is set (optional)
%             default 0.99
%
% OUTPUT nse    : noise estimate (in grayvalues)
%        thresh : threshold for any type of low evidence suppression
%        grad   : buffer with the gradients (allows one to display the 
%                 approximately Rayleigh distributed gradient field)
%
% Matlab-function by KQ 2003 (replaces C-version by GD)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if(nargin == 1)
   c = 0.01;
end;

if(~(isa(img,'double') | isa(img,'uint8')))
   error('Invalid data type entered');
end;

if(isa(img,'uint8'))
   auxImg = double(img)/255;
   [nse,thresh,grad]=NoiseEstim(auxImg,c);
   nse = 255*nse;
   thresh = 255*nse;
   grad = uint8(round(grad*255));
else
   [nse,thresh,grad]=NoiseEstim(img,1-c);
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LOCAL FUNCTION NoiseEstim
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [nse,thresh,grad]=NoiseEstim(auxImg,c)


% find number of dimensions of image
imSize = size(auxImg);
if length(find(imSize-1)) == 1
    nDims = 1;
    auxImg = auxImg(:);
else
    nDims = length(imSize);
end

% filter image. Do gradient, then take only values from 2:end-1, because
% at the edges, forward differences are calculated.
% gradMag is the norm of the gradient vectors.

% the gradient is clearly slower than using conv,conv2 or possibly convn!
centralDifference = [-0.5,0,0.5];
switch nDims
    % code does not work for 1D - I can't seem to get the right results!
%     case 1
%         fx = conv(auxImg,centralDifference);
%         gradMag = abs(fx(2:end-1));
%         
    case 2
        fx = conv2(auxImg,centralDifference','same');
        fy = conv2(auxImg,centralDifference,'same');
        [fx,fy] = gradient(auxImg);
        % calculate gradient magnitude; remove image border
        gradMag = sqrt(fx(2:end-1,2:end-1).^2+fy(2:end-1,2:end-1).^2);
        
    case 3
        
        fx = convn(auxImg,centralDifference','same');
        fy = convn(auxImg,centralDifference,'same');
        fz = convn(auxImg,reshape(centralDifference,[1,1,3]),'same');
        
        % calculate gradient magnitude; remove image border
        gradMag = sqrt(fx(2:end-1,2:end-1).^2+...
            fy(2:end-1,2:end-1).^2+...
            fz(2:end-1,2:end-1).^2);
 
   otherwise
        error('imNoiseEstime needs a 2-3D array as input image')
end

% reclaim memory
[fx,fy, fz, auxImg] = deal([]);

% sort and remove non-zero or NaN elements from gradMag
gradMag = gradMag(:); % make vector
gradMag = sort(gradMag( (gradMag~=0 & ~isnan(gradMag)) ));

% Gradient magnitude is (approximately) Raleigh-distributed. Look 
% for first mode of the distribution in cumulative histogram. 
% To find the mode in, say, 20 data points, take values from 1:10, and
% subtract them from the values from 3:12. The peak of the histogram will
% be where the values lie closest together. Taking the mean of the
% gradients corresponding to the pair will give us an estimate of the peak
% of the histogram. Rinse, repeat for n times.

n = 10;
windowStart = [0,floor([.1:.4/(n-2):.5]*length(gradMag))];
idxVector = [1:floor(length(gradMag)/2)]';

% loop (for memory reasons - it's actually not really slower!)
maxVal = zeros(n-1,1);
for i=1:n-1
    differenceVector = gradMag(idxVector + windowStart(i+1)) - ...
        gradMag(idxVector);
    [minDifference, minIdx] = min(differenceVector);
    % maximum is in the middle between the two bounding values
    maxVal(i) = mean(gradMag(minIdx + windowStart([1,i+1])));
end

% average the maxVals for an estimate of the maximum of the mode
mode = mean(maxVal);

% according to the paper: noise = mode / sqrt(2);
% BUT: since the gradient has been computed over the distance 2 instead of 1
% (lowpass filtering) which is gives a "noise gradient" of a factor 2 too 
% small (the probability to have intensity variation between two locations in an 
% unstructured, noisy image is independent of the distance between the 2 samples
% used to compute the gradient). thus, instead of division by sqrt(2) one has 
% to multiply the mode with 2/sqrt(2)=sqrt(2).
	
switch nDims
    case 2
        multFact = sqrt(2);
    case 3
        % it seems that there is no need to multiply for 3D
        multFact = 1;
end
nse = multFact*mode;
thresh = sqrt(-2.0*log(c))*mode;
grad = gradMag;

% OLD CODE, REPLACED 2/05 by jonas 
% 
% % Convolve the image with the central difference kernel to 
% % obtain the image gradient; calculate gradient magnitude    
% diff=[-0.5;0;0.5];               
% ax = conv2(auxImg, diff', 'same');  
% ay = conv2(auxImg, diff, 'same'); 
% mag = sqrt((ax.*ax) + (ay.*ay)); 
% mag = mag(2:end-1,2:end-1);
% 
% % sort gradient magnitudes in ascending order, 
% % reduce to non-zero elements
% mag = sort(mag(:));
% magr = mag(min(find(mag)):end);
% 
% 
% % Gradient magnitude is (approximately) Raleigh-distributed. Look 
% % for first mode of the distribution in cumulative histogram. 
% winSize = floor([.1:.05:.5]*length(magr));
% 
% diffMat = [1:floor(length(magr)/2)]';
% for i = 1:length(winSize)
%     diffMat = [diffMat, diffMat(:,1)+winSize(i)]; % 10 cols
% end
% data = magr(diffMat);% 10 cols
% diffData = []; 
% for i = 2:length(winSize)+1
%     diffData = [diffData, data(:,i)-data(:,1)]; % 9 cols
% end
% idx = min(diffData); %???
% maxlist = [];
% for i = 1:length(winSize)
%     row = find(diffData(:,i)==idx(i));
%     maxVal = (data(row,i+1)+data(row,1))/2;
%     maxlist = [maxlist;maxVal];
%     % middle of the window position where difference between both borders is minimal
% end
% mode = mean(maxlist);





