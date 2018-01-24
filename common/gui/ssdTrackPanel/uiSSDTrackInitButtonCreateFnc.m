function uiSSDTrackInitButtonCreateFnc
% create function for the init button in the SSD Track panel

panelInfo = get(gcbf,'UserData');

if(~isReadyFg & isempty(panelInfo.stack))
   set(gcbo,'Enable','off');
else
   set(gcbo,'Enable','on');
end;
