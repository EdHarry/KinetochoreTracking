function figH = runUiSSDTrackPanel(viewH)
%RUNUISSDTRACKPANEL starts and runs the SSD tracking panel
%
% SYNOPSIS :  figH = runUiSSDTrackPanel(viewH)
%
% INPUT  viewH : (optional) handle of a view panel which
%                 is connected to the SSD tracking panel
% 
% OUTPUT  figH : handle to the SSD tracking panel

% check whether there exists already a panel
auxH = findobj(0,'Type','figure','Tag','UISSDTRACKPANEL');

if(isempty(auxH))
	% start the panel
   figH = uiSSDTrackPanel;
else
   if(length(auxH) ~= 1)
      error('fatal: more than 1 SSD track panels open');
   else
      figH = auxH(1);
      figure(figH);
   end;
end;

if( nargin == 1)
   % connect to the requested view panel
   panelInfo = get(figH,'UserData');
   panelInfo.connectedView = viewH;
   set(figH,'UserData',panelInfo);
end;
