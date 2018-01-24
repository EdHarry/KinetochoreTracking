function uiSSDTrackGoButtonPressFnc
% callback for the GO button in the SSD Track panel

pI = get(gcbf,'UserData');
success = 1;


if(isempty(pI.stack))
   % presently, the go on! button without loaded stack has the 
   % function of one single match
   uiSSDTrackNextButtonPressFnc(pI,gcbf);
else
   while((~isempty(pI.nextIndx)) & success)
      success = uiSSDTrackNextButtonPressFnc(pI,gcbf);
      pI = get(gcbf,'UserData');
      pause(1);
   end;
end;

   
