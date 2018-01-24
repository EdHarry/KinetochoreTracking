function gbT = LoadFile(gbT)
% callback for the acquisition panel

displ('File will be opened using settings defined in View Panel (Tiff Options)');

% Reads Tiff Option from View Panel
[adepth,algorithm,res,msg]=readTiffOptions;

% Load

gbT.subFrame = [];
gbT.movie.stackIndx = [];
gbT.movie.data = [];
gbT.movie.cMap = [];

fileNameFilter = strcat(gbT.acqDir,filesep,'*.tif');
[fName,dirName] = uigetfile(fileNameFilter,'Acquisition Panel: Load ...');
if( isa(fName,'char') & isa(dirName,'char'))
   % Reads the file
   [gbT.data,adepth,msg,res] = imreadnd([dirName,fName],adepth,algorithm);
else 
   return;
end;

% display the data in ...
gbT.viewPanelH = queryUiViewPanel(gbT.viewPanelH);
if(~isempty(gbT.viewPanelH))
   % ... the connected view panel
   uiViewPanelShowImg(gbT.data,0,gbT.viewPanelH);
else
   % ... in an old existing or a new one
   gbT.viewPanelH = ...
      uiViewPanelShowImg(gbT.data);
end;

[gbT.acqDir,fbody,fno,fext] = getfilenamebody([dirName,fName]);
if(isempty(fno))
   gbT.acqFile = strcat(fbody,'1',fext);
else
   gbT.acqFile = strcat(fbody,int2str(str2num(fno)+1),fext);
end;

% Modifies Tiff options and save menu
modifyTiffOptions(adepth,res,msg);

