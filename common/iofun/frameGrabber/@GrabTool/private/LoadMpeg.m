function gbT = LoadMpeg(gbT)
% callback for the acquisition panel

% get panel information
gbT.subFrame = [];
gbT.movie.stackIndx = [];
gbT.data = [];

% get the default directory
dirEditH = findobj(gcbf,'Tag','ACQ_DIR');

fileNameFilter = strcat(gbT.acqDir,filesep,'*.mpg');
[fName,dirName] = uigetfile(fileNameFilter,'Acquisition Panel: Load MPEG ...');
if( isa(fName,'char') & isa(dirName,'char'))
   [gbT.movie.data, gbT.movie.cMap] = mpgread([dirName,fName]);
else 
   return;
end;

% display the movie in ...
gbT.viewPanelH = queryUiViewPanel(gbT.viewPanelH);
if(isempty(gbT.viewPanelH))
   gbT.viewPanelH = uiViewPanel;
end;
figure(gbT.viewPanelH);
colormap(gbT.movie.cMap);
movie(gca,gbT.movie.data,gbT.movie.repeat,gbT.movie.fps);

gbT.acqDir = getfilenamebody([dirName,fName]);

