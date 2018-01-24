function uiSSDTrackGoButtonCreateFnc
% create function for the GO ON! button in the SSD Track panel

panelInfo = get(gcbf,'UserData');

if(isempty(panelInfo.tImg))
   set(gcbo,'Enable','off');
else
   if(~isempty(panelInfo.stack))
      set(gcbo,'Enable','on');
   end;
end;
