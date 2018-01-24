function uiAcqMenuFileFnc
% callback for the SSD track panel menu

% check if there is a stack loaded
panelInfo = get(gcbf,'UserData');
delH = findobj(gcbo,'Tag','UISSDTRACKMENU_FILE_DELSTACK');

% set the delete menu handle 
if(isempty(panelInfo.stack))
   set(delH,'Enable','off');
else
   set(delH,'Enable','on');
end;
