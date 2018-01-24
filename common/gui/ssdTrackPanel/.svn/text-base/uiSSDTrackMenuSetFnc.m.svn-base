function uiSSDTrackMenuSetFnc
% call back for the paremater set menu in the SSD Track panel

pI = get(gcbf,'UserData');


listStr = {'Full','Shift only','Shift + Scale','Congruent','Shift + Rotation',...
      'Full / linear','Shift only / linear','Shift + Scale / linear',...
      'Congruent / linear','Shift + Rotation / linear'};

[indx,status] = listdlg('ListString',listStr,...
   'Name','SSD Track Parameter Set',...
   'InitialValue',pI.trackCmd.paramsetOpt,...
   'SelectionMode','single','ListSize',[100,100]);

if(status == 1)
   pI.trackCmd.paramsetOpt = indx;
end;

set(gcbf,'UserData',pI);
