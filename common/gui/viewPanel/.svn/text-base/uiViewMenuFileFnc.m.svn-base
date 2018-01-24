function uiViewMenuFileFnc
% callback for the view panel menu

% check if there is image data loaded
axesH = get(gcbf,'CurrentAxes');
imgH = findobj(axesH,'Type','image','Parent',axesH);
saveH = findobj(gcbo,'Tag','UIVIEWMENU_FILE_SAVE');

% get the save menu handle 
if(isempty(imgH))
   set(saveH,'Enable','off');
else
   set(saveH,'Enable','on');
end;
   

