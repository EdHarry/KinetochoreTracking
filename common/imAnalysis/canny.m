function [e,ax,ay,as2,thresh] = canny(varargin)
%CANNY Find edges in an intensity image.
%   
% SYNOPSIS: [e,ax,ay,as2,thresh] = canny(I,thresh,sigma,c);
%
%      The Canny method finds edges by looking for local maxima of the
%      gradient of I. The gradient is calculated after smoothing with 
%      a Gaussian filter. The method uses two thresholds, to detect 
%      strong and weak edges, and includes the weak edges in the output 
%      only if they are connected to strong edges. 
%
% INPUT:
%   I      - image matrix, can be of class uint8, uint16, or double.
%          - alternatively: image name in single quotes, e.g.  
%            'kerato.tif'. Possible formats: 'bmp', 'cur', 'hdf',
%            'ico','jpg' or 'jpeg', 'pcx', `png', 'tif' or 'tiff', 
%            'xwd' (see help for the imread function).
%            
%   optional input variables:
%
%   THRESH - is a two-element vector in which the first element is 
%            the low threshold, and the second element is the high   
%            threshold used in low evidence suppression for the  
%            gradient vector (absolute values). If you specify a scalar
%            for THRESH, this value is used as the high threshold. 
%            If you do not specify THRESH, or if THRESH is empty ([]), 
%            CANNY chooses low and high values automatically.
%   SIGMA  - is the standard deviation of the Gaussian filter. The 
%            default SIGMA is 2; the size of the filter is chosen 
%            automatically, based on SIGMA. 
%   c      - is the confidence probability to which the imNoiseEstim 
%            threshold is set. (The function imNoiseEstim uses Raleigh
%            distribution to calculate the lower threshold). 
%            The default c is 0.99
%
% OUTPUT: 
%   E      - binary edge map
%   AX     - x component of the gradient vector
%   AY     - y component of the gradient vector
%   AS2    - smoothed image
%   THRESH - vector containing the lower and upper thresholds used for 
%            low evidence suppression (see THRESH above) 
%   figure - smoothed image + edge
%            figure is normally supressed. To see it, uncommend the last
%            paragraph of the main function.

% Katharina Quirin, 2002
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[a,thresh,sigma,c] = parse_inputs(varargin{:});

% Transform to a double precision intensity image
if isa(a, 'uint8') | isa(a, 'uint16') 
    a = im2double(a);
end

% The output edge map:
m = size(a,1);
n = size(a,2);
e = repmat(logical(0), m, n);

% Magic numbers
GaussianDieOff = .0001;  
PercentOfPixelsNotEdges = .95; % Used for selecting higher threshold

% Design the filters - a gaussian and its derivative

pw = 1:30;  % possible widths of convolution kernel (depend on sigma)
ssq = sigma*sigma;
width = max(find(exp(-(pw.*pw)/(2*sigma*sigma))>GaussianDieOff));
if isempty(width)
    width = 1;  % the user entered a really small sigma
end
t = (-width:width);

gau = exp(-(t.*t)/(2*ssq)).'/(2*pi*ssq);  % the gaussian 1-D filter
                                          % kernel is column vector
x=sum(gau);
gau = gau./x; % normalized gaussian 1-D filter: weights in the kernel 
              % sum to 1, which avoids increasing or decreasing the 
              % average grey-level when the kernel is used for smoothing
diff=[-0.5;0;0.5]; 
              % central difference kernel to calculate the image gradient
              % (normalized)

% Convolve the filters with the image in each direction. The canny edge 
% detector first requires convolutions with the gaussian, and then with 
% the derivitave of a gaussian. Instead of convolving the filters first,  
% I will first smooth the image by convolution with the normalized 
% gaussian filter in order to estimate the noise from the smoothed image 
% Make a call to conv2 to do the convolution  down each column.

as1 = conv2(a, gau, 'same');
   % convolution in y direction, then 
as2= conv2(as1, gau', 'same'); 
   % convolution in x direction, (transpose kernel)

% figure; imshow(a,[]);
% hold on; title ('raw image')
% figure; imshow(as2,[]);
% hold on; title ('smoothed image')

% Estimate background noise using the function NoiseEstim in order to 
% derive the lower threshold for the low evidence suppression 
disp('arrive at imNoiseEstim')
nse = imNoiseEstim(as2,c); 
disp('imNoiseEstim OK')

% Convolve the smoothed image with the central difference kernel to 
% obtain the image gradient     
ax = conv2(as2, diff', 'same');  
ay = conv2(as2, diff, 'same'); 
% recall: diff is a column vector
mag = sqrt((ax.*ax) + (ay.*ay)); 

magmax = max(mag(:));
if magmax>0
    magn = mag / magmax;  
end         % normalized gradient vector will be used only to                 
            % calculate higher threshold for the low evidence  
            % suppression later on


% Select the thresholds                                                                      
if isempty(thresh) 
    [counts,x]=imhist(magn, 64);
    hthresh = min(find(cumsum(counts)>PercentOfPixelsNotEdges*m*n))/64;
    highThresh=hthresh*magmax;    
            % absolute value for minimal gradient vector threshold
    lowThresh = 3*sqrt(2)*nse;    
            % absolute value for maximal gradient vector threshold
    thresh = [lowThresh highThresh];
elseif length(thresh)==1
    highThresh = thresh;
    lowThresh = 3*sqrt(2)*nse   
    thresh = [lowThresh highThresh];
elseif length(thresh)==2
    lowThresh = thresh(1);
    highThresh = thresh(2);
    if (lowThresh > highThresh)%(lowThresh >= highThresh)
        error('Thresh must be [low high]');
    end
end

% The next step is to do the non-maximum supression (NMS) followed 
% by low evidence suppression (LES)
% We will accrue indices which specify ON pixels in strong edgemap
% The array e will become the weak edge map.
 
idxLocalMax = NMS(mag,ax,ay); % non-maximum supression
idxWeak = idxLocalMax(mag(idxLocalMax) > lowThresh); % LES(I)
e(idxWeak)=1;
idxStrong = [idxWeak(mag(idxWeak) > highThresh)];

rstrong = rem(idxStrong-1, m)+1;
cstrong = floor((idxStrong-1)/m)+1;
% calculate row and column coordinates for strong pixels

e = bwselect(e, cstrong, rstrong, 8); % LES(II): 
% use only those weak pixels connected to strong pixels

% 	[i,j] = find(e);
% 	figure; imshow(as2,[]);
% 	hold on; title ('smoothed image + edge');
% 	hold on; h=plot(j,i,'r.');
% 	set(h,'MarkerSize',2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Local Function : NMS
%
function idxLocalMax = NMS(mag,ax,ay)
% non-maximum suppression for matrix mag and gradient vector(ax;ay)

% SYNOPSIS: magmax = nonmax(mag,ax,ay)
%        
% IN:
% mag          - norm of gradient vector
% ax           - x(column) component of gradient vector
% ay           - y(row) component of gradient vector
%
% OUT:
% idxLocalMax  - vector with indices of pixels with local maximum in mag

re = size(mag,1);

% divide ax and ay by norm of gradient vector, so that norm of 
% resulting gradient vector will be 1
toDoIndx = find(mag ~= 0);
ax(toDoIndx) = ax(toDoIndx)./mag(toDoIndx);
ay(toDoIndx) = ay(toDoIndx)./mag(toDoIndx);


% calculate  column (x) and row (y) indices of the points 
% displaced by +1 along the gradient vector:
magpc = (floor(([1:length(mag(:))]-1)/re)+1)'+1*ax(:);
magpr = (rem([1:length(mag(:))]-1, re)+1)'+1*ay(:);

% calculate  column (x) and row (y) indices of the points 
% displaced by -1 along the gradient vector:
magmc = (floor(([1:length(mag(:))]-1)/re)+1)'-1*ax(:);
magmr = (rem([1:length(mag(:))]-1, re)+1)'-1*ay(:);

% prepare table for comparison of mag values:
% c1: mag; c2: calculated mag at distance +1 along the vector gradient
% c3: calculated mag at distance -1 along the vector gradient
Comp=[mag(:), interp2(mag,magpc,magpr,'linear'),...
              interp2(mag,magmc,magmr,'linear')];
% NOTE: if indices of displaced points are out of range, interp2
% will return NaN. therefore, the corresponding (border) pixels will
% never occur in idxLocalMax

% accrue indices of pixels with local maximum in mag:
idxLocalMax = find((Comp(:,1)>Comp(:,2))&(Comp(:,1)>Comp(:,3)));



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Local Function : parse_inputs
%
function [a,thresh,sigma,c] = parse_inputs(varargin)
% OUTPUTS:
%   a      Image Data
%   thresh Threshold value
%   sigma  standard deviation of Gaussian
%   c      confidence probability to which the imNoiseEstim 
%          threshold is set

error(nargchk(1,4,nargin));

I = varargin{1};
if ischar(I),
    a = imread(I);
else a=I;
end

% Defaults
sigma = 2;
thresh = [];
c = 0.99;

% Now get the rest of the arguments 
threshSpecified = 0; % Threshold is not yet specified
sigmaSpecified = 0;  % Sigma is not yet specified
for i = 2:nargin
    if prod(size(varargin{i}))==2 & ~threshSpecified
        thresh = varargin{i};
        threshSpecified = 1;
    elseif isempty(varargin{i}) & ~threshSpecified
        % Thresh = [];
        threshSpecified = 1;   
    elseif prod(size(varargin{i}))==1 
        if  threshSpecified
            if ~sigmaSpecified
                sigma = varargin{i};
                sigmaSpecified = 1; 
            else
                c = varargin{i};
            end
        elseif ~threshSpecified
            thresh = varargin{i};
            threshSpecified = 1;
        else
            error('Invalid input arguments');
        end
    end
end


if sigma<=0
    error('Sigma must be positive');
end

if c<=0
    error('c must be positive');
end

if isrgb(I)
    error('RGB images are not supported. Call RGB2GRAY first.');
end