function gbT = GrabToFile(gbT)
% callback for the grab & save button in the acquisition panel

showCheckH = findobj(gcbf,'Tag','SHOW_ACQ');

% grab first ...
[s,gbT.data] = mexFgInterface('grab');
leftV=squeeze(gbT.data(:,:,1))';

% ... then eventually show ...
if(get(showCheckH,'Value'))
   if(~isempty(gbT.viewPanelH))
      gbT.viewPanelH = ...
         uiViewPanelShowImg(leftV,0,gbT.viewPanelH);
   else
      gbT.viewPanelH = ...
         uiViewPanelShowImg(leftV,0);
   end;
end;

gbT.movie.stackIndx = [];
gbT.movie.data = [];
gbT.movie.cMap = [];

% ... then save to disk
dirEditH = findobj(gcbf,'Tag','ACQ_DIR');
fileEditH = findobj(gcbf,'Tag','ACQ_FILE');
dirName = get(dirEditH,'String');
fileName = get(fileEditH,'String');
if(isempty(fileName))
   errordlg('no filename specified','Acquisition Panel: Error');
   return;
end;
[filePath, fileBody, fileNo, fileExt] = getFilenameBody(fileName);
if(isempty(fileBody) | isempty(fileExt) )
   errordlg('invalid filename specified','Acquisition Panel: Error');
   return;
end;

fullName = strcat(dirName,'\',fileName);
imwrite(leftV,fullName,'Compression','none');

% update the file edit control
if(isempty(fileNo))
   newFileName = strcat(fileBody,'1',fileExt);
else
   newFileName = strcat(fileBody,int2str(str2num(fileNo)+1),fileExt);
end;
set(fileEditH,'String',newFileName);
