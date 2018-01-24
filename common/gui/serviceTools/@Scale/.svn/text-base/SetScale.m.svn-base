function SetScale(sc)
%SCALE\SETSCALE setup a scale
%
% SYNOPSIS SetScale(scale)
%
% INPUT  scale: an scale object
%
% c: 12/8/99	dT

% Check if there's an option window open
auxH=findobj(0,'Type','figure','Tag','SCALE_OPTIONS');
if (isempty(auxH))
   auxH=OptionWindow(sc);
   % Connect the Scale to the options window   
   set(auxH,'UserData',sc);
else
	figure(auxH);
end;
