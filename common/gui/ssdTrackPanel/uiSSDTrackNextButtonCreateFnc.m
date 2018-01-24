function uiSSDTrackNextButtonCreateFnc
% create function for the next button in the SSD Track panel

panelInfo = get(gcbf,'UserData');

if(isempty(panelInfo.tImg))
   set(gcbo,'Enable','off');
else
   set(gcbo,'Enable','on');
end;
