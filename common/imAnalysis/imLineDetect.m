function [resp,ori,nse,scaleMap,maxMap] = imLineDetect(img,scales,linetype,conf)
%IMLINEDETECT line filter / detector according to 
%  Koller et al. ICCV'95
%  Danuser PhD.  1997
%  (unpublished modifications 1998)
%
%  SYNOPSIS [resp,ori,nse,scaleMap,maxMap] = imLineDetect(img,scales,linetype,conf)
%
%  INPUT img : image structure
%              *.data 2-d image matrix  (double or uint8 supported)
%              *.perm permutation status 'C' or 'M'
%        scales: scale vector on which the multi-scale filtering is performed
%        linetype: definition of the linetypes sought
%    	             is constructed as bitwise ored of the following types
%                   1 positive line (bright line on dark background)
%                   2 negative line (dark line on bright background)
%                   4 wave line (line with bright and dark stripe on 
%                                grey background)
%                   therefore:
%                   3 looks for positive AND negative lines
%                   5 looks for positive AND wave lines
%                   6 looks for negative AND wave lines
%                   7 looks for ALL three types
%       conf: (optional) confidence level on which the thresholds are set 
%             (default 0.99)
%
%  OUTPUT resp : filter response  (double map)
%         ori  : line orientation [0,pi]; perpendicular to line direction
%         nse  : noise level on which data thresholding is performed
%         scaleMap : contains the indexes of the chosen scales
%                    best visualised by
%                    imshow(scalemap,codeColors(max(max(scaleMap))));
%         maxMap : (optional) contains a pixel coded with 1 for all the local maxima
%                  0 otherwise

% check the data type image data
if(isempty(img.data))
   error('no image data available');
end;

if(isa(img.data,'uint8'))
   img.data = double(img.data)/255;
end;
if(img.perm == 'M')
   img.data = permute(img.data,[2,1]);
   img.perm = 'C';
end;

% build the options structure
if(length(scales)>0)
   opt.scales = scales;
else
   error('invalid scale vector');
end;
if((linetype > 0) & (linetype < 8))
   opt.linetype = linetype;
else
   error('invalid linetype specified');
end;

if(nargin < 4)
   conf = 0.99;
end;

[aux,opt.noise] = imNoiseEstim(img.data,1-conf);

if(nargout == 5)
   opt.nonMaxSupp = 1 ;
else
   opt.nonMaxSupp = 0 ;
end;

% start line filtering / detection 
[resp,ori,dtls] = mexLineDetect(img.data,opt);

% permute response and orientation fields
resp = permute(resp,[2,1]);
ori = permute(ori,[2,1]);

% fill other return fields
nse = dtls.noise;
scaleMap = permute(dtls.scaleSelect,[2,1]);

if(~isempty(dtls.maxMap))
   maxMap = permute(dtls.maxMap,[2,1]);
end;

return;