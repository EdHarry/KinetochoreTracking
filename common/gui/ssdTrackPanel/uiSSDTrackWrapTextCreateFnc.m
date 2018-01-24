function uiSSDTrackWrapTextCreateFnc
% create function for the wrap checkbox text in the SSD Track panel

panelInfo = get(gcbf,'UserData');

if(isempty(panelInfo.tImg))
   set(gcbo,'Enable','off');
else
   set(gcbo,'Enable','on');
end;
