function uiViewMenuAnalyzeFnc
% callback for the view panel menu

% check if there is image data loaded
axesH = get(gcbf,'CurrentAxes');
imgH = findobj(axesH,'Type','image','Parent',axesH);
subMenuH = get(gcbo,'Children');

% set the work menu handle 
if(isempty(imgH))
   for i=1:length(subMenuH)
      set(subMenuH,'Enable','off');
   end;
else
   for i=1:length(subMenuH)
      set(subMenuH,'Enable','on');
   end;
end;
