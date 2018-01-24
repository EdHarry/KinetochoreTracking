function figH = uiViewPanelShowImg(img, new, exH)
%UIVIEWPANELSHOWIMG opens (if necessary) and shows an image in the 
% view panel
% 
% SYNOPSIS uiViewPanelShowImg(img, new, exH)
% 
% INPUT 	img : image (must be 2-dimensional matrix)
%        new : (optional) opens a new view panel independent of the
%              existence of another (default: 0 -> search for an
%              existing view panel)
%        exH : (optional) existing handle; used instead of 
%               opening or searching a uiViewPanel
%
% OUTPUT figH : figure handle (if [] something went wrong)
%

% define minimal window size which should be pertained
minDim = uiViewPanelGetDflt('minWinDim');
titleHeight = uiViewPanelGetDflt('winTitleHeight');

% check input image
if( isempty(img) | ~((ndims(img) == 2) | (ndims(img) == 3)))
   error('no image entered');
end;
imgDim = size(img);

% get the figure handle
switch(nargin)
case 1,
   figH = findUiViewPanel;
   if (isempty(figH))
      figH=uiViewPanel;
   end;
case 2,
   figH = findUiViewPanel;
   if (isempty(figH) | new == 1)
      figH=uiViewPanel;
   end;
case 3,
   figH = findUiViewPanel(new,exH);
   if (isempty(figH))
      figH=uiViewPanel;
   end;
end;
figH = figure(figH);

% set dimensions of the axes
axesH = get(figH,'CurrentAxes');
axesPos = get(axesH,'Position');
axesPos(3) = imgDim(2);
axesPos(4) = imgDim(1);
set(axesH,'Position',axesPos);

% check whether window is large enough
figPos = get(figH,'Position');
nomDim(1) = axesPos(3) + 20;
nomDim(2) = axesPos(4) + 45;
if(nomDim(1) < minDim(1))
   nomDim(1) = minDim(1);
end;
if(nomDim(2) < minDim(2));
   nomDim(2) = minDim(2);
end;

% make sure that the entire window is on the screen
screenDim = get(get(figH,'Parent'),'ScreenSize');
if((figPos(1) + nomDim(1)) > screenDim(3))
   figPos(1) = screenDim(3) - nomDim(1);
end;
if((figPos(2) + nomDim(2) + titleHeight) > screenDim(4))
	figPos(2) = screenDim(4) - nomDim(2) - titleHeight;
end;

figPos = [figPos(1:2),nomDim];
set(figH,'Position',figPos);

% the gamma correction tool requires double type images. 
% If the input image is uint8 then convert it here for display 
% purposes
if(isa(img,'uint8'))
   img = double(img)/255;
else
    img = double(img)/255;
    minI = min(min(img));
    maxI = max(max(img));
    
    if( ~((minI >= 0) & (maxI <= 1)) )
        % the incoming data is NOT normalized and thus we have to take care of it
        if( (minI < 0) )
            % the assumption is that some mathematical data is shown; thus stretch it to the 
            % min/max range
            img = (img - minI) / (maxI - minI);
        else
            % the assumption is that some image data of an unknown number of bits (8 bit at minimum) 
            % is displayed
            bitRange = ceil(log2(maxI));
            if bitRange < 8
                % enforce a minimum of 8 bit
                bitRange = 8;
            end;
            img = img / 2^bitRange;
        end;    
    end;
end;

if(get(0,'ScreenDepth') == 8)
   global nColors__ ;
   imshow(img,nColors__);
else
   imshow(img);
end;

% apply gamma correction if a gamma correction window is associated with it
imT = ImAdjustTool(figH);
ctrlWinH = QueryCtrlWindowHandle(imT);
if(~isempty(ctrlWinH))
   cmap = QueryCMap(imT);
   if(~isempty(cmap))
      colormap(cmap);
   end;
end;

VPanelResize(figH,[0 20 0 0]);
% image(img);
% set(axes,'UserData',img);