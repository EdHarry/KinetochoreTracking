function gbT= OpenGrabPanel(gbT)
%GRABTOOL/OPENGRABPANEL opens a new image acquisition window
%
% SYNOPSIS gbT = OpenGrabPanel(gbT)
%
% INPUT  gbT: an grabTool object
%
% c: 24/8/99	dT

% 
if (nargout ~= 1)
   error('fetch the grabtool object');
end;

% Check if there's a corresponding grab window open
auxH = findobj(0,'Tag','GRABWINDOW');

if (isempty(auxH))
   [gbT.fg.typeList,gbT.fg.nr]=InitFrameGrabbers(gbT);
   % Open the grab window window
   gbT.grabPanelH=GrabPanel(gbT);
   % Connect the gbTool to the grab window   
   set(gbT.grabPanelH,'UserData',gbT);
else
   figure(auxH);
end;