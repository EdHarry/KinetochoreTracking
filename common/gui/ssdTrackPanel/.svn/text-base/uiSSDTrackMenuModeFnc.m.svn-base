function uiSSDTrackMenuModeFnc
% call back for the paremater mode menu in the SSD Track panel

pI = get(gcbf,'UserData');

listStr = {'No Statistics','Statistics','Diagnostics'};

[indx,status] = listdlg('ListString',listStr,...
   'Name','SSD Track Computation Mode',...
   'InitialValue',pI.trackCmd.compOpt+1,...
   'SelectionMode','single','ListSize',[100,100]);

if(status == 1)
   pI.trackCmd.compOpt = indx-1;
end;

%% At the moment the diagnostics SW is not reviewd; thus 
%% set the index back to statistics until this is accomplished
%
% the disganostics option should now work Sept-17-1998
%if(pI.trackCmd.compOpt == 2)
%   pI.trackCmd.compOpt = 1;
%end;
      
set(gcbf,'UserData',pI);
