function gbT = LoadClippedStack(gbT)

displ('Enter the minimal and maximal intensity in bit (or just press Okay for a simple load)');

% Dialog
prompt={'gmin (bit):','gmax (bit):'};
def={'',''};
dlgTitle='Define stack global boudaries';
lineNo=1;
answer=inputdlg(prompt,dlgTitle,lineNo,def);

if isempty(answer)
   displ('Action interrupted');
   return
end

% Conversion
gmin=2^str2num(char(answer(1)))-1;
gmax=2^str2num(char(answer(2)))-1;

% Check for entered parameters
if isempty(gmin) | isempty(gmax)
   % Function switch
   ft=1;
else
   % Display
   displ(strcat('Range=[',num2str(gmin),':',num2str(gmax),']'));
   % Function switch
   ft=0;
end

oldZoomStatus = 0;
oldPixelStatus = 0;

gbT.movie.data = [];
gbT.movie.cMap = [];

% get default directory
fileNameFilter = strcat(gbT.acqDir,filesep,'*.tif');
[fName,dirName] = uigetfile(fileNameFilter,...
   'Acquisition Panel: Load Stack ...');

% quit function if filename is not valid
if( ~(isa(fName,'char') & isa(dirName,'char')))
   return;
end;

% else read first image
auxI = imread([dirName,fName]);
typ=class(auxI);

% display the data in ...
gbT.viewPanelH = queryUiViewPanel(gbT.viewPanelH);
if(~isempty(gbT.viewPanelH))
   % ... the connected view 
   gbT.viewPanelH = ...
      uiViewPanelShowImg(auxI,0,gbT.viewPanelH);
else
   % ... in an old existing or a new one
   gbT.viewPanelH = uiViewPanelShowImg(auxI);
end;

% make the connected view panel the active figure
figH = figure(gbT.viewPanelH);

% toggle the zoom and pixel functions if necessary
zoomMenuH = findobj(figH,'Tag','UIVIEWMENU_TOOLS_ZOOM');
if(strcmp(get(zoomMenuH,'Checked'),'on'))
   % toggle the zoom menu status to off
   uiViewMenuZoomFnc(zoomMenuH,figH);
   oldZoomStatus = 1;
else
   pixMenuH = findobj(figH,'Tag','UIVIEWMENU_TOOLS_PIXEL');
   if(strcmp(get(pixMenuH,'Checked'),'on'))
      % toggle the pixel menu status to off
      uiViewMenuPixelFnc(pixMenuH,figH);
      oldPixelStatus = 1;
   end;
end;

% set the text in the view panel
textH = findobj(figH,'Type','uicontrol','Tag','UIVIEWTEXT');
oldText = get(textH,'String');
oldTextStatus = get(textH,'Visible');
set(textH,'String','crop region with left button');
set(textH,'Visible','on');

% interactively crop the image
[auxI,gbT.subFrame] = imcrop;

% set the old text back
set(textH,'Visible',oldTextStatus);
set(textH,'String',oldText);

% restore the old view panel status
if(oldPixelStatus)
   % toggle the pixel menu status back to on
   uiViewMenuPixelFnc(pixMenuH,figH);
else
   if(oldZoomStatus)
      % toggle the zoom menu status back to on
      uiViewMenuPixelFnc(zoomMenuH,figH);
   end;
end

% read the stack
switch ft
   case 0
[gbT.data,msg,res,gbT.movie.stackIndx,nInSeries,gbT.subFrame] = ...
   imreadstacknd([dirName,fName],gmin,gmax,gbT.movie.maxImgInStack, gbT.subFrame );
case 1
[gbT.data,gbT.movie.stackIndx,nInSeries,gbT.subFrame] = ...
   imreadstack([dirName,fName],gbT.movie.maxImgInStack, gbT.subFrame );
      % Conversion of the loaded data to double for manipulation
      %switch class(gbT.data)
      switch typ
      case 'uint8'
         gbT.data=double(gbT.data)/(2^8-1);
      case 'uint16'
         gbT.data=double(gbT.data)/(2^16-1);
      otherwise
         error('Bit depth not supported');
      end
      % Output
      msg='Stack loaded without normalization';
otherwise
   error('Error encountered');
end

nImg = length(gbT.movie.stackIndx);

% display the data in ...
gbT.viewPanelH = queryUiViewPanel(gbT.viewPanelH);
if(~isempty(gbT.viewPanelH))
   % ... the connected view panel
   uiViewPanelShowImg(gbT.data(:,:,1),0,...
      gbT.viewPanelH);
else
   % ... in an old existing or a new one
   gbT.viewPanelH = ...
      uiViewPanelShowImg(gbT.data(:,:,1));
end;

% set the new filename
[gbT.acqDir,fbody,fno,fext] = getfilenamebody([dirName,fName]);
gbT.acqFile = strcat(fbody,int2str(gbT.movie.stackIndx(nImg)+1),fext);

% Display result
displ(msg);
