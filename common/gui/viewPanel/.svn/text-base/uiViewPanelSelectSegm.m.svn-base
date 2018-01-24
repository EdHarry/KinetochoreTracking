function [x,y] = uiViewPanelSelectSegm(exH)
%UIVIEWPANELSELECTSEGM service function for the View Panel 
% allowing one to interactively select a segment
% 
% SYNOPSIS [x,y] = uiViewPanelSelectSegm(exH)
%
% INPUT  
%          exH : external handle (optional)
%
% OUTPUT   x,y : coordinate of the segment end points 

x = [];
y = [];

if(nargin == 1)
   ownH = queryUiViewPanel(exH);
   if(isempty(ownH))
      return;
   end;
else
   ownH = gcbf;
end;

oldZoomStatus = 0;
oldPixelStatus = 0;

% toggle the zoom and pixel functions in view panel if necessary
zoomMenuH = findobj(ownH,'Tag','UIVIEWMENU_TOOLS_ZOOM');
if(strcmp(get(zoomMenuH,'Checked'),'on'))
   % toggle the zoom menu status to off
   uiViewMenuZoomFnc(zoomMenuH,ownH);
   oldZoomStatus = 1;
else
   pixMenuH = findobj(ownH,'Tag','UIVIEWMENU_TOOLS_PIXEL');
   if(strcmp(get(pixMenuH,'Checked'),'on'))
      % toggle the pixel menu status to off
      uiViewMenuPixelFnc(pixMenuH,ownH);
      oldPixelStatus = 1;
   end;
end;

% set the text in the view panel
textH = findobj(ownH,'Type','uicontrol','Tag','UIVIEWTEXT');
oldText = get(textH,'String');
oldTextStatus = get(textH,'Visible');
set(textH,'String','select segment with left button');
set(textH,'Visible','on');

% select segment interactively;
[x,y] = getSegment(ownH);

% set the old text back
set(textH,'String',oldText);
set(textH,'Visible',oldTextStatus);

% restore the old view panel status
if(oldPixelStatus)
   % toggle the pixel menu status back to on
   uiViewMenuPixelFnc(pixMenuH,ownH);
else
   if(oldZoomStatus)
      % toggle the zoom menu status back to on
      uiViewMenuZoomFnc(zoomMenuH,ownH);
   end;
end;
