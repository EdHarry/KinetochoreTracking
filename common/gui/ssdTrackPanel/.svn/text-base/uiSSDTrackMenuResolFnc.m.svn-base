function uiSSDTrackMenuResolFnc
% call back for the resolution menu in the SSD Track panel

pI = get(gcbf,'UserData');

prompt = {'Shift resolution [Pixel]','Shape resolution multiplicator'};
title = 'SSD Track computational resolution';
lineNo = 1;
strVal1 = sprintf('%f',pI.trackCmd.resol(1));
strVal2 = sprintf('%f',pI.trackCmd.resol(2));
dfltVals = {strVal1,strVal2};
result = inputdlg(prompt,title,lineNo,dfltVals);

if(isempty(result))
   return;
else
   pI.trackCmd.resol(1) = str2num(result{1});
   pI.trackCmd.resol(2) = str2num(result{2});
   set(gcbf,'UserData',pI);
end;