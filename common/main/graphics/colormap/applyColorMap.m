function [img3C,pixClasses]=applyColorMap(img,scores,range,cmap,convFactor)
% applyColorMap applies a colormap to score classes in a user-defined range.
%
% This function expects the range to be either
%   i)  symmetric around 0 (e.g. [-15 15])
%   ii) completely positive (e.g. [0 15]).
% The scores must have the same distribution.
%
% SYINOPSIS   [img3C,pixClasses]=applyColorMap(img,scores,range,cmap,convFactor)
%
% INPUT       img        : image to which the color-coded scores will be
%                          overlaid.
%             scores     : matrix of scores. It must have the same dimensions
%                          of img.
%             range      : score range: [lowerBound upperBound]. 
%                          Accepted ranges: [-bound bound] | [0 bound].
%                          Pass range=[] to use the default range [-15 15].
%             cmap       : colormap to be used. Pass cmap=[] to use the
%                          default colormap. The default colormap is meant
%                          to represent the range [-15 15], where the
%                          0-class corresponds to the 16th entry in the
%                          colormap. cmap can be either a colormap matrix
%                          (nx3 matrix), or the name of a colormap (string).
%                          Type 'help graph3d' to see a number of useful colormaps.
%             convFactor : (optional) This factor stretches the matrix
%                          scores to belong to the passed range. convFactor is 
%                          by default 1.
% 
% REMARK      The spanned range and the number of colors in the colormap
%             must match.
%             For instance, the range [-15 15] spans 31 classes, which is
%             the number of colors in the default colormap (it is a 31x3
%             matrix)
%
% OUTPUT      img3C      : color-coded RGB image
%             pixClasses : matrix containing (at each position) a pointer
%                          to the corresponding colormap index.
%
% Aaron Ponti, November 25th, 2002

% Default input parameter values
dfltRange=[-15 15];
dfltConvFactor=1;
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
             1.0000    1.0000    1.0000; % White -> No color
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
if nargin<4 | nargin>5
    error('The function expects 4 or 5 parameters.');
end
if nargin==4 | (nargin==5 & isempty(convFactor))
    convFactor=dfltConvFactor;
end

% Check range
range=sort(range); % Make sure that lowerBound<upperBound
if isempty(range)
    range=dfltRange;
    classCase=1;
else
    if range(1)==0 & range(2)>0
        classCase=2;
    elseif sign(range(1))~=sign(range(2)) & abs(range(1))==abs(range(2))
        classCase=1;  
    else
        fprintf(1,'Range not valid. Using the default range [%d %d].',dfltRange);
        range=dfltRange;
        classCase=1;
    end
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
    
% Check that range and colormap match
if length([range(1):range(2)])~=size(cmap,1)
    error('Range and colormap dimensions do not match');
end
    
% Render img lighter (if needed)
minImgI = min(img(:));
maxImgI = max(img(:));
if minImgI==0.5 & maxImgI==1
    % The image was already treated
elseif minImgI < maxImgI
    img=0.5+nrm(img,1)/2;  
end

% Make sure the scores are in the passed range and are 'integer'
scores=round(scores*convFactor);
scores(find(scores<range(1)))=range(1);
scores(find(scores>range(2)))=range(2);

switch classCase
    case 1 % Range symmetric around 0
        % Assign scores to classes
        pixClasses=scores+1+fix(size(cmap,1)/2);
        % Find positions where scores do not belong to the 0-class (e.g. 16 for the range [-15 15])
        [y x]=find(pixClasses~=fix(size(cmap,1)/2)+1);
    case 2
        % Assign scores to classes (no shift)
        pixClasses=scores;
        % Find positions where scores do not belong to the 0-class (1)
        [y x]=find(pixClasses~=0);
    otherwise
        error('This case is not supported');
end


% Create 3 channels
img3C(:,:,1)=img; img3C(:,:,2)=img; img3C(:,:,3)=img;

% Apply color
for k=1:length(y)
    img3C(y(k),x(k),1)=img3C(y(k),x(k),1)*cmap(pixClasses(y(k),x(k)),1);
    img3C(y(k),x(k),2)=img3C(y(k),x(k),2)*cmap(pixClasses(y(k),x(k)),2);
    img3C(y(k),x(k),3)=img3C(y(k),x(k),3)*cmap(pixClasses(y(k),x(k)),3);
end