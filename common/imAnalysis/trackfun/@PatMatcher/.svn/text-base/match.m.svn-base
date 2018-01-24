function [pm,pos,info]=Match(pm,tPos,pos0,pSze,cmd)
%MATCH matches a pattern in a search image by minimizing an SSD.
%
%SYNOPSIS : [pm,pos,status]=Match(pm,tPos,pos0,pSze,cmd)
%
% INPUT 
%       tPos: position of the patch center 
%       pos0: first guess of the position in the search image
%       pSze: patch size [width, height, depth] 
%      
%       cmd: command structure with the fields
%            *.mask  : uint8 binary mask to select valid pixels in template(same size as patch)
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
%            *.prec  : computational resolution [res_shift, mult_res_shape]
%            *.compOpt : computational mode
%                        0 = no statistics; 1 = statistics; 
%                        2 = diagnostics;
%            *.interpolOpt : interpolation mode
%                            1 = nearest neighbour; 2 = bilinear; 3 = cubic; 4 =spline
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
% OUTPUT pm : PatMatcher object
%        pos   : position of the matched pattern
%        status


pm.cmd=cmd;
pm.pSze = pSze;
pm.tPos = tPos;
pm.pos0 = pos0;
pos = NaN;
pm.model.params=[];
pm.info.status='Ok';

% add one to patch size for gradient computation
hSze=floor(pm.pSze/2)+1;
pm.patCenter=hSze+1;

% Check if patch (p center) is chosen correctly
if(~all(clipnegative(pm.tPos-hSze)))
   pm.info.status='patch size too large';
   info=pm.info;
   return;
end;

% prepare n-dimensional template and interpolation
pm.dims=ndims(pm.tempImg.data);
crop='crImg(';
intp='interpn(pm.srchImg.data,';
for dim = 1:pm.dims,
   sdim=num2str(dim);
   crop=strcat(crop,['(crCent(' sdim ')-hSze(' sdim ')):(crCent(' sdim ')+hSze(' sdim ')),']);
   intp=strcat(intp,['pm.model.coord(' sdim ',:),']);
end;
intp=strcat(intp,'pm.interp{pm.cmd.interpol})');
crop=crop(1:(length(crop)-1));
crop=strcat(crop,')');

% Cut template from template image
crImg=pm.tempImg.data;
crCent=round(pm.tPos);
pm.template.values = eval(crop);
pm.template.mean=mean(pm.template.values(:));
pm.template.std=std(pm.template.values(:));
pm.template.size=size(pm.template.values);

% Create coordinates for patch and apply mask if not empty
tSze=pm.template.size;
totTSze=length(pm.template.values(:));

if(isempty(pm.cmd.mask))
	[tempCoord{1:pm.dims}]=ind2sub(pm.template.size,1:totTSze);
else
   if all(size(pm.cmd.mask)==pm.template.size)
      pm.template.values=pm.template.values.*pm.cmd.mask;
      [tempCoord{1:pm.dims}]=ind2sub(pm.template.size,find(pm.cmd.mask));
   else
      pm.info.status='bad mask';
      return;
	end;
end;

for l = 1:pm.dims
	pm.template.coord(l,:)=tempCoord{l};
end;

pm.patch=zeros(tSze);
% Set correct crop size
hSze=hSze-1;
crCent=pm.patCenter;
% coord # of center pixel
cent=ceil(prod(tSze)/2);

% Set intial parameter change
dpar=[];
%Counter for iteration
iterCt=0;
while(isempty(dpar) | (any(abs(dpar(:,iterCt).*osz_mask)>pm.cmd.prec) & iterCt<pm.cmd.maxIter) )
   iterCt=iterCt+1;
   %Compute new coords according to model
   pm=model(pm);
   %add offset from search image
	for dim = 1:pm.dims
      pm.model.coord(dim,:)=pm.model.coord(dim,:)+pm.pos0(dim);
   end;
   %Compute interpolated values and store in patch
   pm.patch(:)=eval(intp);
   % Check if all coords are in searchimg
   if(any(isnan(pm.patch)))
      pm.info.status='patch out of search image'; 
      info=pm.info;
      return;
   end;
   %Intensity correction
   pm=Radiometric(pm);
   
   % Check if coords out of excursion
   if(norm(pm.template.coord(:,cent)+pm.pos0'-pm.patCenter'-pm.model.coord(:,cent))>pm.cmd.maxEx)
      pm.info.status='patch exceeded max excursion';
      info=pm.info;
      return;
   end;

   %Comp Grad(I_pat)
   [dI{1:pm.dims}]=gradient(pm.patch);
   dI=[dI(2) dI(1) dI(3:pm.dims)];
   dIp=[];
   %Ignore the frame voxels
   for l = 1:pm.dims
      crImg=dI{l};
      crI=eval(crop);
      dIp=[dIp crI(:)];
   end;
   
   crImg=pm.template.values;
   It=eval(crop);
   crImg=pm.patch;
   Ip=eval(crop);

   % error with current params: e = mA*params-mB
   mB=It-Ip;
	mA=[];
   for l=1:size(dIp,1)
      mA=[mA ; dIp(l,:)*pm.model.grad(:,:,l)];
   end;
	dpar(:,iterCt)= mA\mB(:);
   %Check for oscillation
   if(iterCt>3)
	   %osz_mask=abs(dpar(:,iterCt))<abs(dpar(:,iterCt-1));
	   osz_mask=abs(dpar(:,iterCt))<abs(dpar(:,iterCt-1)) | abs(dpar(:,iterCt-1))<abs(dpar(:,iterCt-2));
   else
      osz_mask=ones(length(dpar),1);
   end;
   pm.model.params=pm.model.params+(dpar(:,iterCt).*osz_mask);
end;
if (iterCt>=pm.cmd.maxIter)
   pm.info.status='maximal iteration depth exceeded';
end;
info = pm.info;
% find center pixel:
pos=pm.model.coord(:,cent);