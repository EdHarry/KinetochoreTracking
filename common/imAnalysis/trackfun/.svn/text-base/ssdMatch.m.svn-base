function [pos,dtls,tImgNew,sImgNew]=ssdMatch(tImg,sImg,tPos,pos0,shape0,pSze,cmd)
%SSDMATCH matches a pattern in a search image by minimizing an SSD.
% In its current version, the function is wrapper for the 
% dLsm library.
%
%SYNOPSIS : [pos,dtls,tImgNew,sImgNew]=ssdMatch(tImg,sImg,tPos,pos0,shape0,pSze,cmd)
%
% INPUT tImg: template image structure
%       sImg: search image structure
%             Both images have the same structure with the fields
%             *.data : 2d image matrix (uint8 and double supported)
%             *.perm : permutation status
%                      'C' or 'M'
%             *.prefilter : 1 prefilter the image
%                         : 0 leave as it is
%       tPos: position of the template center 
%       pos0: first guess of the position in the search image
%       shape0: first guess of the affine transformation matrix [mx0, sx0; sy0, my0]
%               (optional field: can also be empty)
%       pSze: patch size [width, height] -> will be converted
%             into # of pixels
%       cmd: command structure with the fields
%            *.mask  : uint8 binary mask to select valid pixels in template
%            *.lori  : direction enforced on displacement estimate (obsolete if paramsetOpt is set 
%                      to an unconstrained geometric model). 
%                      The value is defined in a left handed image coordinate system,
%                      as the angle between the x1-axis and the direction perpendicular to the linear
%                      pattern.
%                      If lori == [] and a linear constraint
%                      is applied to the displacement estimation lori is calculated from 
%                      analyzing the texture in the template.
%            *.sigma : parameter for prefiltering
%            *.maxEx : maximal excursion in position or shape
%            *.resol : computational resolution [res_shift, mult_res_shape]
%            *.compOpt : computational mode
%                        0 = no statistics; 1 = statistics; 
%                        2 = diagnostics;
%            *.interpolOpt : interpolation mode
%                            1 = nearest neighbour; 2 = bilinear;
%            *.paramsetOpt : parameter set
%                            1 = all parameters
%                            2 = shift only
%                            3 = shift + scales
%                            4 = congruent
%                            5 = shift + rotation
%                            6 = all parameters with line displacement constraint
%                            7 = shift only with line displacement constraint
%                            8 = shift + scales with line displacement constraint
%                            9 = congruent with line displacement constraint
%                           10 = shift + rotation with line displacement constraint
%            *.patchSeriesOpt : store patch series during iteration
%                               1 = yes; 0 = no;
%            *.verbose : run ssd matcher in verbose mode
%                        1 = yes; 0 = no;
%
% OUTPUT pos : position of the matched pattern
%        dtls : structure containing details of the match with the 
%               following fields:
%               *.status: status report from the dLsm library
%               *.lori: local orientation enforced on displacement estimate
%                       The value is defined in a left handed image coordinate system,
%                       as the angle between the x1-axis and the direction perpendicular to the linear
%                       pattern.
%               *.shape: shape matrix of the estimated affine transform
%               *.goffit: goodness of fit evaluation by [sigma0, cross_correl_coefficient]
%               *.preci: precision of the estimated parameters
%               *.nPix: number of pixels used for tracking
%               *.patchStack: stack with patches acquired during iteration
%             
%        tImgNew : updated template image structure
%        sImgNew : updated search image structure
%                  both structures are equal to the input images;
%                  they are always of type "double" and the 
%                  permutation flag is set to 'C'. Thus, for 
%                  visualization the images must be permuted back.

global dfltsFilename__;

% check if there is data 
if(isempty(tImg.data) | isempty(sImg.data))
   if(isempty(tImg.data))
      dtls.status = -10; % match status with dLsmErrcds.h file
   else
      dtls.status = -11; % match status with dLsmErrcds.h file      
   end;
   dtls.shape = [1,0;0,1];
   dtls.lori = [];
   pos = pos0;
   tImgNew.data = [];
   tImgNew.perm = 'M';
   tImgNew.prefilt = 1;
   sImgNew.data = [];
   sImgNew.perm = 'M';
   sImgNew.prefilt = 1;
   return;
end;

% check the data type of the images
if(isa(tImg.data,'uint8'))
   tImg.data = double(tImg.data)/255;
end;
if(isa(sImg.data,'uint8'))
   sImg.data = double(sImg.data)/255;
end;
if(~isempty(cmd.mask))
   if(isa(cmd.mask,'uint8'))
      cmd.mask = permute(cmd.mask,[2,1]);
   else
      cmd.mask = [];
   end;   
end;

% permute if necessary
if(tImg.perm == 'M')
   aux = permute(tImg.data,[2,1]);
   tImg.data = aux;
   tImg.perm = 'C';
end;
if(sImg.perm == 'M')
   aux = permute(sImg.data,[2,1]);
   sImg.data = aux;
   sImg.perm = 'C';
end;

% building the command structure for the ssdHandler
ssdCmd.fname = [];  % change that later
ssdCmd.tImg = tImg.data;
ssdCmd.tMsk = cmd.mask;
ssdCmd.sImg = sImg.data;
ssdCmd.pSze = pSze;
ssdCmd.tPos = tPos;
ssdCmd.pos0 = pos0;
ssdCmd.shape0 = shape0;
ssdCmd.maxEx = cmd.maxEx;
ssdCmd.resol = cmd.resol;
ssdCmd.cParams = cmd.lori;
ssdCmd.opts(1) = cmd.compOpt;
ssdCmd.opts(2) = cmd.interpolOpt;
ssdCmd.opts(3) = cmd.paramsetOpt;
ssdCmd.opts(4) = tImg.prefilter;
ssdCmd.opts(5) = sImg.prefilter;
ssdCmd.opts(6) = cmd.patchSeriesOpt;
ssdCmd.opts(7) = cmd.verbose;

if(tImg.prefilter | sImg.prefilter)
   ssdCmd.filtParam = cmd.sigma;
end;

% call ssd matcher
[status, result] = mexSSDHandler(ssdCmd);

% fill the return values
pos = result.pos;
dtls.lori = result.lori;
dtls.status = result.status;
dtls.shape = result.shape;
dtls.preci = result.preci;
for(i=1:length(dtls.preci))
   if(dtls.preci(i) == 0)
      dtls.preci(i) = NaN;
   end;
end;
dtls.goffit = result.goffit;

if(~isempty(ssdCmd.tMsk))
   dtls.nPix = sum(sum(ssdCmd.tMsk));
else
   dtls.nPix = (pSze(1)+1)*(pSze(1)+2);
end;

if(~isempty(result.patchStack))
   dtls.patchStack = permute(result.patchStack,[2,1,3]);
else
   dtls.patchStack = [];
end;

% get the return images ready
tImgNew.data = result.tImg;
tImgNew.perm = 'C';
tImgNew.prefilter = 0;
sImgNew.data = result.sImg;
sImgNew.perm = 'C';
sImgNew.prefilter = 0;

return;