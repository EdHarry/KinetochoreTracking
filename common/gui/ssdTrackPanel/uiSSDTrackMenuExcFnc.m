function uiSSDTrackMenuExcFnc
% call back for the excursion menu in the SSD Track panel

pI = get(gcbf,'UserData');

prompt = {'Maximal excursion [Pixel]'};
title = 'SSD Track maximal excursion';
lineNo = 1;
strVal1 = sprintf('%d',pI.trackCmd.maxEx);
dfltVals = {strVal1};
result = inputdlg(prompt,title,lineNo,dfltVals);

if(isempty(result))
   return;
else
   pI.trackCmd.maxEx = str2num(result{1});
   set(gcbf,'UserData',pI);
end;

%% pI.trackCmd