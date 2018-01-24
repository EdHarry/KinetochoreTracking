function [newTpos,posRef,maskParams,mask]=uiSSDTrackPanelSetTmask(img,tPos,tDim,maskType,fH)
% service function which sets a mask over the template upon request

posRef = [0,0];

switch maskType
case 1, [posRef] = manSlctCtr(img,tPos,tDim);
   newTpos = tPos;
   mask = [];
   maskParams = [];
case 2, [newTpos,maskParams,mask]=manSlctSegm(img,tPos,tDim,fH);
case 3, [newTpos,maskParams]=getDICLine(img,tPos,tDim);
   mask = [];
case 4, [newTpos,maskParams,mask]=getDICLineAndMask(img,tPos,tDim);
case 5, [newTpos,maskParams]=getPosLine(img,tPos,tDim);
   mask = [];	
case 6, [newTpos,maskParams,mask]=getPosLineAndMask(img,tPos,tDim);
otherwise,
   mask=[];
   maskParams=[];
   newTpos = tPos;
end;

%---------------------------------------------------------------------------------
%
% Local functions for mask computation
%
%---------------------------------------------------------------------------------

function [posRef] = manSlctCtr(img,tPos,roiDim);


figH = figure;
set(figH,'Name','SSD Track Template');
cRect = [round(tPos)- roiDim/2,roiDim];
tImg = imcrop(img,cRect);
imshow(tImg);

msgbox('Select reference point','SSD Track request','modal');

pos = getPt(figH);
posRef = cRect(1:2)-1+pos-tPos;
close(figH);

%---------------------------------------------------------------------------------
function [pos,alpha,mask]=manSlctSegm(img,tPos,tDim,figH)
global dicLineWidth__;


% plot the selected box 
figure(figH),
hold on;
plotrect([tPos-tDim/2,tDim],'y-');
hold off;

boundsX = [];
while(length(boundsX)<2)
   msgbox('Select segment','SSD Track request','modal');
   [segmX, segmY] = uiViewPanelSelectSegm(figH);
   [boundsX,boundsY] = segmentIntersectRect(segmX,segmY,[tPos-tDim/2,tDim]);
end;

% replot the image to erase the temporary template border
uiViewPanelShowImg(img,0,figH);

% compute new center position 
pos = [mean(boundsX), mean(boundsY)];

% compute orientation of the line (direction normal to it)
alpha = atan2(boundsY(2)-boundsY(1), boundsX(2)-boundsX(1)) + pi/2;
if(alpha > pi) 
   alpha = alpha - pi;
end;
if(alpha < 0)
   alpha = alpha + pi;
end;

% compute the new mask
% 1st: get the segment bounds relative to the new template position
[boundsX,boundsY] = segmentIntersectRect(segmX,segmY,[pos-tDim/2,tDim]);

mask = uint8(zeros(tDim(2)+1,tDim(1)+1));
mSze = size(mask);
mCtr = tDim/2+1;

% fill the mask
cCtr = round(pos);
for i = -ceil(dicLineWidth__/2):1:ceil(dicLineWidth__/2)
   xS = [boundsX(1),boundsY(1)] - cCtr + i*[cos(alpha),sin(alpha)];
   xE = [boundsX(2),boundsY(2)] - cCtr + i*[cos(alpha),sin(alpha)];
   pts = bresenham(xS,xE);
   for(j=1:size(pts,2))
      mPos = (pts(:,j))' + mCtr;
      if((mPos(1)>0)&(mPos(1)<=mSze(2))&(mPos(2)>0)&(mPos(2)<=mSze(1)))
         mask(mPos(2),mPos(1)) = uint8(1);
      end;
   end;
end;


%---------------------------------------------------------------------------------

function [pos,alpha]=getDICLine(img,pos0,roiDim)
global dicLineScales__;
global lineConfProb__;

cDim = roiDim + 2*ceil(4*dicLineScales__(end));
cRect = [round(pos0)-cDim/2,cDim];
I.data = imcrop(img,cRect);
ISze = size(I.data);
if( ((ISze(1) - 1)~=cRect(4)) | ((ISze(2) - 1)~=cRect(3)))
   error('invalid image cropping');
end;
I.perm = 'M';

[resp,ori,nse] = imLineDetect(I,dicLineScales__,4,lineConfProb__);
[d,alpha] = lineHough(resp,ori);

pos = cRect(1:2) + cDim/2 + [cos(alpha),sin(alpha)]*d;

%---------------------------------------------------------------------------------

function [pos,alpha,mask]=getDICLineAndMask(img,pos0,roiDim)
global dicLineScales__;
global lineConfProb__;
global debuggingMode__;

d = 100;   % dummy initial value for d
nIter = 0;

cDim = roiDim + 2*ceil(4*dicLineScales__(end));
cCtr = cDim/2+1;

% place the cropped image at the right place
while ((abs(d) > 1) & (nIter < 3))
   cRect = [round(pos0)-cDim/2,cDim];
   I.data = imcrop(img,cRect);
   ISze = size(I.data);
   if( ((ISze(1) - 1)~=cRect(4)) | ((ISze(2) - 1)~=cRect(3)))
      error('invalid image cropping');
   end;
   I.perm = 'M';
   
   [resp,ori,nse,scales] = imLineDetect(I,dicLineScales__,4,lineConfProb__);
   [d,alpha,hMap] = lineHough(resp,ori);
   
   if((nIter == 0) & (debuggingMode__ == 1))
      figure;
      imshow(-1*hMap,[]);
   end;
   
   pos0 = cRect(1:2) + cDim/2 + [cos(alpha),sin(alpha)]*d;
   nIter = nIter + 1;
end;

pos = pos0;

mask = uint8(zeros(roiDim(2)+1,roiDim(1)+1));
mSze = size(mask);
mCtr = roiDim/2+1;
% compute the central axis' end points
alphaC = mod(alpha+pi/2,pi);
alpha0 = atan(mSze(1)/mSze(2));
if((alphaC>alpha0) & (alphaC<= pi - alpha0))
   auxX(2) = (mSze(1)+1)/2;
   if(alphaC == pi/2)
      auxX(1) = 0;
   else
      auxX(1) = auxX(2)/tan(alphaC);
   end;
else
   auxX(1) = (mSze(2)+1)/2;
   auxX(2) = auxX(1)*tan(alphaC);
end;
% interpolate the scales on the central line
nSmpl = round(6*sqrt(sum(auxX.^2))) + 1;   % sampling rate 0.33
[sclX,sclY,sclV] = improfile(scales,[cCtr(1)-auxX(1),cCtr(1)+auxX(1)],...
   [cCtr(2)-auxX(2),cCtr(2)+auxX(2)],nSmpl,'nearest');

% fill the mask
for i = 1:length(sclX)
   if(sclV(i))
      % width is given as 1/0.25 * sigma_opt
      xS = [sclX(i),sclY(i)]-cCtr+2*dicLineScales__(sclV(i))*[cos(alpha),sin(alpha)];
      xE = [sclX(i),sclY(i)]-cCtr-2*dicLineScales__(sclV(i))*[cos(alpha),sin(alpha)];
      pts = bresenham(xS,xE);
      for(j=1:size(pts,2))
         mPos = (pts(:,j))' + mCtr;
         if((mPos(1)>0)&(mPos(1)<=mSze(2))&(mPos(2)>0)&(mPos(2)<=mSze(1)))
            mask(mPos(2),mPos(1)) = uint8(1);
         end;
      end;
   end;
end;

if(debuggingMode__ == 1)
   fh = gcf;
   figure;
   imshow(resp,[]);
   figure;
   imshow(ori,[]);
   figure(fh);
end;
   
%---------------------------------------------------------------------------------

function [pos,alpha]=getPosLine(img,pos0,roiDim)
global posLineScales__;
global lineConfProb__;

cDim = roiDim + 2*ceil(4*posLineScales__(end));
cRect = [round(pos0)-cDim/2,cDim];
I.data = imcrop(img,cRect);
ISze = size(I.data);
if( ((ISze(1) - 1)~=cRect(4)) | ((ISze(2) - 1)~=cRect(3)))
   error('invalid image cropping');
end;
I.perm = 'M';

[resp,ori,nse] = imLineDetect(I,posLineScales__,1,lineConfProb__);
[d,alpha] = lineHough(resp,ori);

pos = cRect(1:2) + cDim/2 + [cos(alpha),sin(alpha)]*d;

%---------------------------------------------------------------------------------

function [pos,alpha,mask]=getPosLineAndMask(img,pos0,roiDim)
global posLineScales__;
global lineConfProb__;
global debuggingMode__;

d = 100;   % dummy initial value for d
nIter = 0;

cDim = roiDim + 2*ceil(4*posLineScales__(end));
cCtr = cDim/2+1;

% place the cropped image at the right place
while ((d > 1) & (nIter < 3))
   cRect = [round(pos0)-cDim/2,cDim];
   I.data = imcrop(img,cRect);
   ISze = size(I.data);
   if( ((ISze(1) - 1)~=cRect(4)) | ((ISze(2) - 1)~=cRect(3)))
      error('invalid image cropping');
   end;
   I.perm = 'M';
   
   [resp,ori,nse,scales] = imLineDetect(I,posLineScales__,1,lineConfProb__);
   [d,alpha] = lineHough(resp,ori);
   
   pos0 = cRect(1:2) + cDim/2 + [cos(alpha),sin(alpha)]*d;
   nIter = nIter + 1;
end;

pos = pos0;

mask = uint8(zeros(roiDim(2)+1,roiDim(1)+1));
mSze = size(mask);
mCtr = roiDim/2+1;
% compute the central axis' end points
alphaC = mod(alpha+pi/2,pi);
alpha0 = atan(mSze(1)/mSze(2));
if((alphaC>alpha0) & (alphaC<= pi - alpha0))
   auxX(2) = (mSze(1)+1)/2;
   if(alphaC == pi/2)
      auxX(1) = 0;
   else
      auxX(1) = auxX(2)/tan(alphaC);
   end;
else
   auxX(1) = (mSze(2)+1)/2;
   auxX(2) = auxX(1)*tan(alphaC);
end;
% interpolate the scales on the central line
nSmpl = round(6*sqrt(sum(auxX.^2))) + 1;   % sampling rate 0.33
[sclX,sclY,sclV] = improfile(scales,[cCtr(1)-auxX(1),cCtr(1)+auxX(1)],...
   [cCtr(2)-auxX(2),cCtr(2)+auxX(2)],nSmpl,'nearest');

% fill the mask
for i = 1:length(sclX)
   if(sclV(i))
      % width is given as 1/0.417 * sigma_opt
      xS = [sclX(i),sclY(i)]-cCtr+1.2*posLineScales__(sclV(i))*[cos(alpha),sin(alpha)];
      xE = [sclX(i),sclY(i)]-cCtr-1.2*posLineScales__(sclV(i))*[cos(alpha),sin(alpha)];
      pts = bresenham(xS,xE);
      for(j=1:size(pts,2))
         mPos = (pts(:,j))' + mCtr;
         if((mPos(1)>0)&(mPos(1)<=mSze(2))&(mPos(2)>0)&(mPos(2)<=mSze(1)))
            mask(mPos(2),mPos(1)) = uint8(1);
         end;
      end;
   end;
end;

if(debuggingMode__ == 1)
   fh = gcf;
   figure;
   imshow(resp,[]);
   figure;
   imshow(mask,[]);
   figure(fh);
end;
