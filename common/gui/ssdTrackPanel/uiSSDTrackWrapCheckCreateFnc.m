function uiSSDTrackWrapCheckCreateFnc
% create function for the wrap checkbox in the SSD Track panel

panelInfo = get(gcbf,'UserData');

if(isempty(panelInfo.tImg))
   set(gcbo,'Enable','off');
else
   set(gcbo,'Enable','on');
end;

set(gcbo,'Value',0);
