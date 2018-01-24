function uiSSDTrackLogEditButtonPressFnc
% callback for the file browser button

% get handles of the edit field
fileEditH = findobj(gcbf,'Tag','UISSDTRACKLOGEDIT');
[fpath,fbody,fno,fext] = getfilenamebody(get(fileEditH,'String'));
initFileString = strcat(fpath,filesep,'*.log');

% get the new filename
[newFileString,newDirString] = ...
   uiputfile(initFileString,'SSD Track Panel: Set Logfile');

if( isa(newFileString,'char') & isa(newDirString,'char') )
   % set the new data	
   set(fileEditH,'String',[newDirString,newFileString]);
end;



