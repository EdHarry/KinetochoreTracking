function uiSSDTrackMenuLoadStackFnc
% callback for the Load Stack menu in the SSD track panel

% get panel information
pI = get(gcbf,'UserData');
pI.tImg = [];
pI.sImg = [];
pI.tPos = [];
pI.tPosOff = [];
pI.tDispOri = [];
pI.tPosRef = [0,0];
pI.pos0 = [];
pI.pos = [];
pI.posOff = [];
pI.patchDim = [];
pI.stackIndx = [];
pI.stackTot = [];
pI.nImgGrabbed = 0;

% get default directory and filename
fileEditH = findobj(gcbf,'Tag','UISSDTRACKLOGEDIT');
[fpath,fbody,fno,fext] = getfilenamebody(get(fileEditH,'String'));
dfltFileName = strcat(fpath,filesep,'*.tif');
[fName,dirName] = uigetfile(dfltFileName,...
   'SSD Track Panel: Load Stack ...');
if( isa(fName,'char') & isa(dirName,'char'))
   pI.stackName = [dirName,fName];
   [pI.stack,pI.stackIndx,pI.stackTot] = ...
      imreadstack(pI.stackName,...
      uiSSDTrackPanelGetDflt('maxImgInStack'));
else 
   return;
end;

% set the new filename
fName = uiSSDTrackPanelGetDflt('logfname');
set(fileEditH,'String',[dirName,fName]);

%% the following block, immediately starting template 
%% initialization is commented out. It is more sensical to 
%% force the user to press the INIT button

if(length(pI.stackIndx)>1)
   pI.nextIndx = 2;
   % start the template init function
   
   %% this is obsolete in the new version of the SSDtrackPanel 
   %% as template initialization is performed only by pressing the 
   %% the INIT button
   %  [pI.tImg,pI.tPos,pI.tPosOff,pI.tPosRef,pI.patchDim,pI.trackCmd,pI.connectedView]=...
   %      uiSSDTrackPanelInitTimg(pI.stack(:,:,1),pI.stackIndx(1),...
   %      pI.tImgMaskType,...
   %      pI.tImg,...
   %      pI.tPos,...
   %      pI.tPosRef,...
   %      pI.patchDim,...
   %   	pI.trackCmd,...
   %      pI.connectedView);
   % set the displacement orientation in the template coordinate frame 
   % pI.tDispOri = pI.trackCmd.lori;
   %   pI.pos0 = pI.tPos;
   %   pI.posOff = pI.tPosOff;
   % initPropagator('ID');
   % setPropagatorPosEstimate(pI.tPos+pI.tPosOff-1);
else
   pI.stack = [];
end;


if(~isempty(pI.stack))
   
   % show the first image
   pI.connectedView=uiViewPanelShowImg(pI.stack(:,:,1),1,pI.connectedView);

   %%  the two following lines have been moved to 
   %%  uiSSDTrackInitButtonPressFnc()
   %%  as in the new version template initialization does no more 
   %%  take place automatically after loading a stack.
   % set(findobj(gcbf,'Tag','UISSDTRACKWRAPCHECK'),'Enable','on');
   % set(findobj(gcbf,'Tag','UISSDTRACKWRAPTEXT'),'Enable','on');
   set(findobj(gcbf,'Tag','UISSDTRACKINITBUTTON'),'Enable','on');
   set(findobj(gcbf,'Tag','UISSDTRACKFROMLOGCHECK'),'Enable','on');

   % set the text of the view panel
   textH = findobj(pI.connectedView,'Type','uicontrol','Tag','UIVIEWTEXT');
   set(textH,'String',sprintf('%d/%d',pI.stackIndx(1),...
      pI.stackIndx(1)-1+pI.stackTot));
   set(textH,'Visible','on');
end;

% assign all the data to the figure
set(gcbf,'UserData',pI);

