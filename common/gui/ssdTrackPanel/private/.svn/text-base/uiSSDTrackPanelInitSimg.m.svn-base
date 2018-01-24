function [sImg,newPos0,newPosOff] = ...
   uiSSDTrackPanelInitSimg(img,pos0,posOff,pDim,trackCmd)
% service function for the initialization of the search image

ctrPos = posOff + pos0 - 1;

cDim = pDim + 2*trackCmd.maxEx + ...
      2*ceil(uiSSDTrackPanelGetDflt('bdFactor')*trackCmd.sigma);	

cRect = [round(ctrPos) - cDim/2,cDim];

if(isempty(img))
   sImg.data = grab(cRect);
else
   % crop the image
   sImg.data = imcrop(img,cRect);
end;
sImgSze = size(sImg.data);
if( ((sImgSze(1) - 1)~=cRect(4)) | ((sImgSze(2) - 1)~=cRect(3)))
   sImg.data = [];
end;
sImg.perm = 'M';
sImg.prefilter = 1;
newPosOff = cRect(1:2);
newPos0 = ctrPos - newPosOff + 1;