function [binMask,resp,ori] = imLineExtract(img,opt)
%IMLINEEXTRACT line extraction calling imLineDetect(), 
% which is a wrapper for the old C - line filter implemented by Danuser. 
% The extractor function applies nonmax suppression on the filter output
% and hysteresis thresholding according to Canny, PAMI 1986
% in order to extract a pixelated line.
%
% SYNOPSIS [binMask,resp,ori] = imLineExtract(img,opt)
%
% INPUT img : double or uint8 matrix with greyvalue image data
%       opt : (optional) data structure with options for expert users 
%             if the whole structure or some specific fields are 
%             missing, the function sets default values
%
%         *.scales : scale vector for multi scale edge extraction 
%                    default 1:3
%         *.linetype : polarity of extracted line. For index definition 
%                      see help page of imLineDetect()
%                      default 1 = bright line on dark background
%         *.relThresh : relative thresholds for hysteresis thresholding in the 
%                       range [0,1]
%                       [min,max] : the effective thresholds are the computed 
%                                   based on minThresh = min% of all non-zero responses 
%                                            maxThresh = max% of all non-zero responses
%                       if min > 1, the lower threshold will be calculated as
%                                            min * noise(responseMap), 
%                                            where the noise estimate is claculated
%                                            using imNoiseEstim()
%                       if only a scalar value in the range [0,1] is set, then the 
%                       the thresholds will be set to 
%                                            maxThresh = scalar% of all non-zero responses
%                                            0.4 * maxThresh
%
%                       default is max% = 0.8; min = 4*noise(responseMap);
%
% OUTPUT binMask : [0,1] map with 1s at extracted line locations (linels)
%                        and 0s elsewhere 
%        resp :    with filter response at extracted linels
%        ori  :    with loacl orientation of the extracted linels 
%                 (perpendicular to line direction)
%
% SEE ALSO imLineDetect, imOverlayMask (for visualization puproses)


% manage input data and fill in default values where necessary
if nargin == 1
   opt.scales = 1:3;
   opt.linetype = 1;
else
   if ~isfield(opt,'scales')
      opt.scales = 1:3;
   end;
   if ~isfield(opt,'linetype')
      opt.linetype = 1;
   end;
end;

% convert img to the data structure required for imLineDetect()
img.data = img;
img.perm = 'M';
[resp,ori,dummy,dummy,maxMap] = imLineDetect(img,opt.scales,opt.linetype,1);
binMask = repmat(logical(uint8(0)),size(img.data,1),size(img.data,2));


% start hysteresis thresholding

% first step is to set the thresholds
imgSize = size(img.data,1)*size(img.data,2);

% magic numbers 
maxPrct = 0.8;
thresholdRatio = 0.4;
nTimesNseLevel = 4;

srtResp=sort(resp(resp>0));

if ~isfield(opt,'relThresh') 
   highThresh = srtResp(round(maxPrct*length(srtResp)));
   % estimate the noise in the response field to determine the lower threshold
   nseResp = imNoiseEstim(resp);
   lowThresh = nTimesNseLevel*nseResp;
elseif length(opt.relThresh)==1
   if opt.relThresh>=1
      error('The threshold must be less than 1.');
   end
   highThresh = srtResp(round(opt.relThresh*length(srtResp)));
   lowThresh = thresholdRatio*highThresh;
elseif length(opt.relThresh)==2
    highThresh = srtResp(round(opt.relThresh(2)*length(srtResp)));
    if opt.relThresh(1) < 1;
       lowThresh = srtResp(round(opt.relThresh(1)*length(srtResp))); 
    else
      % lower threshold is set based on the noise in the response field
      nseResp = imNoiseEstim(resp);
      lowThresh = opt.relThresh(1)*nseResp;
   end;
   if (lowThresh >= highThresh) | (opt.relThresh(2) >= 1)
      error('lower threshold gets larger the high threshold, or higher threshold > 1.');
   end
end

% second: select two index arrays, one for the weak and one for the strong linels
idxLocMax = find(maxMap);
idxWeak = idxLocMax(resp(idxLocMax) > lowThresh);
idxStrong = idxLocMax(resp(idxLocMax) > highThresh);
% binMask contains all values that are above lowThresh
binMask(idxWeak) = 1; 

% call a morphological operator that accepts all objects (in this case line segments) 
% that touch a strong linel
rstrong = rem(idxStrong-1, size(resp,1))+1;
cstrong = floor((idxStrong-1)/size(resp,1))+1;
binMask = bwselect(binMask, cstrong, rstrong, 8);

% reduce the response and orientation maps
idxBinMask = find(binMask);
auxMap = zeros(size(resp));
auxMap(idxBinMask) = resp(idxBinMask);
resp = auxMap;
auxMap = repmat(NaN,size(img.data,1),size(img.data,2));
auxMap(idxBinMask) = ori(idxBinMask);
ori = auxMap;



