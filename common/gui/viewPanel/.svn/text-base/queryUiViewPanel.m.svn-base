function vpH = queryUiViewPanel(exH)
% service function which queries a handle to a view panel
%
% INPUT exH : handle to presumably existing view panel
% 
% OUTPUT vpH : handle to the figure.
%

vpH = [];

if(isempty(exH))
   return
end;

% check whether the optional handle is a view panel
auxH = findobj(0,'Type','figure','Tag','UIVIEWPANEL');
for( i = 1:length(auxH))
   if(auxH(i) == exH)
         vpH = exH;
   end;
end;
