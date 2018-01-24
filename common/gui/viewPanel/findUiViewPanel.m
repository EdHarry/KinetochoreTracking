function figH = findUiViewPanel(new,exH)
%FINDUIVIEWPANEL service call to find an existing view panel 
% on the system or, if requested, generates a new one.
%
% INPUT new : (optional) generate new panel (default : 0)
%       exH : (optional) external handle to presumably existing
%             view panel
% 
% OUTPUT figH : handle to the figure, [] if nothing found
%

figH = [];

% set new to its default value
if(nargin == 0)
   new = 0;
end;

if(nargin == 2)
   % check whether the optional handle is a view panel
   if(isempty(queryUiViewPanel(exH)))
      figH = [];
   else
      figH = exH;
   end;
end;

if( isempty(figH) )
   % check whether there is an existing view panel
   auxH = findobj(0,'Type','figure','Tag','UIVIEWPANEL');
   if(length(auxH) > 0)
      figH = auxH(1);
   else
      if( new > 0 )
         figH = uiViewPanel;
      else
         figH = [];
      end;
   end;
end;
