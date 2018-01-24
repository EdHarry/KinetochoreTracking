function gbT = LoadFullStack(gbT)
% callback for the Load Stack menu in the acquisition panel

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

% Loads
gbT.subFrame = [];
gbT.movie.data = [];
gbT.movie.cMap = [];

% get default directory
fileNameFilter = strcat(gbT.acqDir,filesep,'*.tif');
[fName,dirName] = uigetfile(fileNameFilter,...
   'Acquisition Panel: Load Stack ...');
if( isa(fName,'char') & isa(dirName,'char'))
   switch ft
   case 0
      displ('Evaluating... please wait');
      [gbT.data,msg,res,gbT.movie.stackIndx] = ...
         imreadstacknd([dirName,fName],gmin,gmax,gbT.movie.maxImgInStack);
   case 1
      auxI = imread([dirName,fName]);
      typ=class(auxI);
      clear auxI;
      [gbT.data,gbT.movie.stackIndx] = ...
         imreadstack([dirName,fName],gbT.movie.maxImgInStack);
      % Conversion of the loaded data to double for manipulation
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
      error('Error encountered.');
   end
else 
   return;
end;
nImg = length(gbT.movie.stackIndx);

% display the data in ...
gbT.viewPanelH = queryUiViewPanel(gbT.viewPanelH );
if(~isempty(gbT.viewPanelH ))
   % ... the connected view panel
   uiViewPanelShowImg(gbT.data(:,:,1),0, gbT.viewPanelH );
else
   % ... in an old existing or a new one
   gbT.viewPanelH  = uiViewPanelShowImg(gbT.data(:,:,1));
end;

% set the new filename
[gbT.acqDir,fbody,fno,fext] = getfilenamebody([dirName,fName]);
gbT.acqFile = strcat(fbody, int2str(gbT.movie.stackIndx(nImg)+1),fext);

% Display result
displ(msg);
