function [img3C,pixClasses]=imDataMapOverlay(img,dataM,range,cmap,opacity)
% imDataMapOverlay Overlay a colored data map to an image.
%
% SYINOPSIS :  
%    [img3C,pixClasses]=imDataMapOverlay(img,dataM)
%    [img3C,pixClasses]=imDataMapOverlay(img,dataM,range,cmap,opacity)
%
% INPUT :
%    img     : image to which the color-coded dataM will be
%              overlaid.
%    dataM   : The data map, a 2D matrix. It must have the same dimensions
%              of img. In an image area where no data is available,
%              the original image is displayed without color being assigned.
%              Set dataM to 'NaN' at these pixels.
%    range   : The range of data (format [lowerBound upperBound]) where
%              colors are assigned according to a specific color
%              map. The color index of 'lowerBound' is the first index while
%              The color index of 'upperBound' is the last index.
%              Out of range data are assigned the min
%              (<lowerBound) and max (>upperBound) color
%              respectively. Default is the whole data range.
%              Pass range=[] to use the default range.
%    cmap    : colormap to be used. Pass cmap=[] to use the
%              default colormap. 
%              Type 'help graph3d' to see a number of useful colormaps.
%    opacity : The level of opacity (0~1). When it is 0, it is
%              fully transparent. When it is 1, no image feature
%              can be seen under dataM color. Default value is
%              0.5. 
% 
%
% OUTPUT :
%    img3C      : Color-coded RGB image
%    pixClasses : Matrix containing (at each pixel) the colormap index for the
%                 corresponding data at that pixel. It has the same size of
%                 'dataM'.
%
% AUTHOR: Lin Ji, Aug. 30, 2005
%
% This function is modified from Araon's 'applyColormap with the following
% added features:
%    1.) In image area where there is no data (indicated by NaN in 'dataM'),
%        the original image will be displayed without assigned color.
%    2.) Add 'opacity' level option so that the overlay can be opaque at
%        different levels. When opacity == 1, one does not see the 
%        underlying image except for no data area.
%    3.) The spanned range does not have to match the number of colors
%        anymore. 

% Default input parameter values
dfltRange   = [];
dfltOpacity = 0.5;
dfltCMap=[   0         0         1.0000; % Blue
             0         0.0645    0.9677;
             0         0.1290    0.9355;
             0         0.1935    0.9032;
             0         0.2581    0.8710;
             0         0.3226    0.8387;
             0         0.3871    0.8065;
             0         0.4516    0.7742;
             0         0.5161    0.7419;
             0         0.5806    0.7097;
             0         0.6452    0.6774;
             0         0.7097    0.6452;
             0         0.7742    0.6129;
             0         0.8387    0.5806;
             0         0.9032    0.5484;
             1.0000    0.8889         0; % Yellow
             1.0000    0.8254         0;
             1.0000    0.7619         0;
             1.0000    0.6984         0;
             1.0000    0.6349         0;
             1.0000    0.5714         0;
             1.0000    0.5079         0;
             1.0000    0.4444         0;
             1.0000    0.3810         0;
             1.0000    0.3175         0;
             1.0000    0.2540         0;
             1.0000    0.1905         0;
             1.0000    0.1270         0;
             1.0000    0.0635         0;
             1.0000         0         0]; % Red

% Check input parameters
if nargin<2
    error('The function expects at least two parameters.');
end
if nargin > 5
    error('The function can not have more than 5 parameters.');
end
if nargin == 2
   range   = [];
   cmap    = [];
   opacity = [];
end
if nargin == 3 
   cmap    = [];
   opacity = [];
end
if nargin == 4 
   opacity = [];
end

% Check range
if isempty(range)
    range=dfltRange;
elseif ~isnumeric(range) | ...
   (isnumeric(range) & length(range) ~= 2)
   error('The range input is not correct. See help applyColormap.');
elseif range(1) > range(2)
   error(['The range is \[lowerBound, upperBound\] where ' ...
      'lowerBound<=upperBound. See help applyColormap.']);
end

% Check cmap
if ischar(cmap);
    cmap=colormap(cmap);
end
if isempty(cmap)
    cmap=dfltCMap;
end
if size(cmap,2)~=3 | size(cmap,3)~=1
    error('The colormap must be an nx3 matrix where n corresponds to the number of classes spanned by range');
end    
    
%Check opacity.
if isempty(opacity)
   opacity = dfltOpacity;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Assign dataM in the selected range to color classes or index of the
% specified color map
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
numColors  = size(cmap,1);
pixClasses = dataM;

% Find positions where dataM are not NaN. 
numInd = find(~isnan(pixClasses));
minS = min(pixClasses(numInd));
maxS = max(pixClasses(numInd));
if isempty(range)
   s1 = minS;
   s2 = maxS;
else
   s1 = range(1);
   s2 = range(2);
end

pixClasses(numInd(find(pixClasses(numInd)<s1)))=s1;
pixClasses(numInd(find(pixClasses(numInd)>s2)))=s2;
if s1==s2
    %Assign the middle color.
    pixClasses(numInd) = numColors/2;
else
    pixClasses(numInd) = (pixClasses(numInd)-s1)/(s2-s1)*(numColors-1)+1;
end

%Make sure pixClasses are integer after scaling to colormap.
pixClasses(numInd) = round(pixClasses(numInd));

% Find image pixels where no score is availble. 
nanInd = find(isnan(pixClasses));
%[y x]=find(~isnan(pixClasses));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Scale image intensity under no score area (nanInd) to 0~1 so that the image
% under these areas is fully displayed without color overlay. Render img under
% real dataM (numInd) according to the specified opacity level.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
img = double(img);
minImgI     = min(img(nanInd));
maxImgI     = max(img(nanInd));
img(nanInd) = (img(nanInd)-minImgI)/(maxImgI-minImgI);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Scale image intensity under score area (numInd) according to opacity level
% so that the underlying image can be seen to some extent.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
minImgI = min(img(numInd));
maxImgI = max(img(numInd));
img(numInd) = (img(numInd)-minImgI)/(maxImgI-minImgI)*(1-opacity) + opacity;

% Create 3 channels
%img3C(:,:,1)=img; img3C(:,:,2)=img; img3C(:,:,3)=img;
img3C = ones(length(img(:)),3);
img3C(:,1)=img(:); img3C(:,2)=img(:); img3C(:,3)=img(:);

% Apply color
img3C(numInd,1) = img3C(numInd,1).*cmap(pixClasses(numInd),1);
img3C(numInd,2) = img3C(numInd,2).*cmap(pixClasses(numInd),2);
img3C(numInd,3) = img3C(numInd,3).*cmap(pixClasses(numInd),3);
img3C = reshape(img3C,[size(img) 3]);

%for k=1:length(y)
%    img3C(y(k),x(k),1)=img3C(y(k),x(k),1)*cmap(pixClasses(y(k),x(k)),1);
%    img3C(y(k),x(k),2)=img3C(y(k),x(k),2)*cmap(pixClasses(y(k),x(k)),2);
%    img3C(y(k),x(k),3)=img3C(y(k),x(k),3)*cmap(pixClasses(y(k),x(k)),3);
%end
