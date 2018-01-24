function GrabPanelCB(gbT,action)
%GRABTOOL\GRABPANELCB callback for grabPanel

cbObj = gcbo;
cbFig =gcbf;
whoCalled = get(cbObj,'Tag');

%Panel specific functions are implemented below, panel independent
%functions in the private directory

%Handle the action off the appropriate caller object.
switch action
   
	case 'Create'
      switch whoCalled
      	case 'GRABWINDOW'
            %init grab window
         case {'ACQ_DIR','ACQ_FILE'}
	      	if(gbT.fg.nr(1)>0)
   	      	set(cbObj,'Enable','on');
      		else
         		set(cbObj,'Enable','inactive');
            end;
         % keep default values for these:   
         case {'FGTYPE','FILE_MENU','FRAME_MENU','LOAD_FILE','LOADSTACK_FULL','LOADSTACK_CLIPPED','LOAD_MPEG'}
      	otherwise   
      		if(gbT.fg.active>1)
         			set(cbObj,'Enable','on');
      		else
         		set(cbObj,'Enable','off');
      	end;
      end;
      
   case 'Delete'
      switch whoCalled
      	case 'GRABWINDOW'
         	% close the framegrabber and free memory
				if(gbT.fg.active~=1)
   				eval(['mexFgInterface(' 39 'close' 39 ')'],'');
				end; 
      end;
   case 'Act'
      switch whoCalled
         case 'ACQ_DIR'
            gbT.acqDir = get(cbObj,'String');
         case 'ACQ_FILE'
            gbT.acqFile = get(cbObj,'String');
         case 'ACQ_TIME'
            gbT.movie.time = str2num(get(cbObj,'String'));
         case 'BROWSE'
            gbT=Browse(gbT);
            SetPath(gbT);
         case 'FGTYPE'
            gbT=OpenFg(gbT);
      	case 'FILE_MENU'
         	FileMenu(gbT);
      	case 'FRAME_MENU'
         	FrameMenu(gbT);
         case 'FRAMES_SEC'
				gbT.movie.fps = str2num(get(cbObj,'String'));
      	case 'FRAME_SELECT_INTERACTIVE'
         	gbT=SelectInteractive(gbT);
      	case 'FRAME_SELECT_PANEL'
         case 'FULL_FRAME'
				%Set max roi to framegrabber (if an fg is active)- ignore errors
				if (gbT.fg.active~=1)
   				eval(['mexFgInterface(' 39 'grab' 39 ',[1,1,2000,2000])'],'');
				end;        
				%reset subFrame
            gbT.subFrame=[];
      	case 'GRAB_TO_FILE'
         	gbT = GrabToFile(gbT);
				StackControls(gbT,'off');
      	case 'GRAB_TO_VIEW'
         	gbT = GrabToView(gbT);
				StackControls(gbT,'off');
      	case 'GRAB_STACK'
         	gbT=GrabToStack(gbT);
				StackControls(gbT,'on');
      	case 'LOAD_FILE'
         	gbT=LoadFile(gbT);
            SetPath(gbT);
				StackControls(gbT,'off');
      	case 'LOADSTACK_FULL'
            gbT=LoadFullStack(gbT);
            SetPath(gbT);
         	StackControls(gbT,'on');
      	case 'LOADSTACK_CLIPPED'
        		gbT=LoadClippedStack(gbT);
            SetPath(gbT);
				StackControls(gbT,'on');
      	case 'LOAD_MPEG'
         	gbT=LoadMpeg(gbT);
            SetPath(gbT);
            EnableMovieCtls(gbT);
         case 'MOVIE_BUTTON'
         	gbT=ProcessMovie(gbT);
            EnableMovieCtls(gbT);
      	case 'MOVIE_REPEAT'         
				gbT.movie.repeat = str2num(get(cbObj,'String'));
      	case 'SELECT_SUBMENU'
         	SelectSubmenu(gbT);
         case 'SAVE_FILE'
            gbT=SaveFile(gbT);
            SetPath(gbT);
         case 'SAVE_MPEG'
            gbT=SaveMpeg(gbT);
            SetPath(gbT);
         case 'SHUTTER_TIME'
            SetShutterTime(gbT);
         case 'STACK_LIST'
            gbT=StackList(gbT);
         case 'STACK_SLIDER'
            gbT=StackSlider(gbT);
   	end;   
end;

%Reconnect the new data to the panel
set(gbT.grabPanelH,'UserData',gbT);

%--------------------------------------------------------

%Implementation of short grabPANEL specific functions:
%FileMenu
%FrameMenu
%SelectSubmenu
%StackControls
%StackSlider
%StackList
%SetPath
%EnableMovieCtls
%Browse
%OpenFg
%SetShutterTime

%--------------------------------------------------------

function FileMenu(gbT)
saveH = findobj(gcbo,'Tag','SAVE_FILE');
savempgH = findobj(gcbo,'Tag','SAVE_MPEG');

% set the save menu handle 
if(isempty(gbT.data))
   set(saveH,'Enable','off');
else
   set(saveH,'Enable','on');
end;

% set the save mpeg menu handle 
if(isempty(gbT.movie.data))
   set(savempgH,'Enable','off');
else
   set(savempgH,'Enable','on');
end;

%--------------------------------------------


function FrameMenu(gbT)
% check if there is image data loaded
varH = findobj(gcbo,'Tag','FULL_FRAME');

% get the save menu handle 
if(isempty(gbT.subFrame))
   set(varH,'Checked','on');
else
   set(varH,'Checked','off');
end;

%--------------------------------------------

function SelectSubmenu(gbT)

varH = findobj(gcbo,'Tag','FRAME_SELECT_INTERACTIVE');

% check the status of the connected view
if(~isempty(gbT.viewPanelH))
   if(~isempty(queryUiViewPanel(gbT.viewPanelH)))
      set(varH,'Enable','on');
   else	
      gbT.viewPanelH= [];
      set(varH,'Enable','off');
   end;
else	
  set(varH,'Enable','off');
end;

%--------------------------------------------

function StackControls(gbT,request)
% service function for enabling/disabling stack controls
figH=gbT.grabPanelH;
nImg = length(gbT.movie.stackIndx);
% the ones with no interference to acquisition functions
%keyboard;
if(nImg~=0)
   h = findobj(figH,'Tag','STACK_SLIDER');
   if(nImg>2)
	% reset values
	set(h,'Min',gbT.movie.stackIndx(1),...
   	'Max',gbT.movie.stackIndx(nImg),...
   	'SliderStep',[1/(nImg-1),1],...
   	'Value',gbT.movie.stackIndx(1));
   set(h,'Enable',request);
	else
   	set(h,'Enable','off');
   end;
	h = findobj(figH,'Tag','STACK_LIST');
	% reset values
	set(h,'String',...
   	int2str(reshape(gbT.movie.stackIndx,nImg,1)), ...
   	'Value',1);
	set(h,'Enable',request);
	h = findobj(figH,'Tag','MOVIE_BUTTON');
	set(h,'Enable',request);
else
   request='off';
   h = findobj(figH,'Tag','STACK_SLIDER');
	% reset values
	set(h,'Enable',request);
	h = findobj(figH,'Tag','STACK_LIST');
	% reset values
	set(h,'Enable',request);
	h = findobj(figH,'Tag','MOVIE_BUTTON');
	set(h,'Enable',request);
end;

% dependent on wether there is movie data update the button string
if(isempty(gbT.movie.data))
   set(h,'String',gbT.movie.buttonStr(1));
else
   set(h,'String',gbT.movie.buttonStr(2));
end;

% UIACQFPSEDIT	 is also used for acquiring stacks
% therefore no change if the famegrabber capability is available
if(gbT.fg.active==1)
   h = findobj(figH,'Tag','FRAMES_SEC');
   set(h,'Enable',request);
end;

% UIACQREPEATEDIT can be enabled only if movie data is available
% and is disabled together with the other stack functionality
if(strcmp(request,'on'))
   if(~isempty(gbT.movie.data))
      h = findobj(figH,'Tag','MOVIE_REPEAT');
      set(h,'Enable',request);
   end;
else
   h = findobj(figH,'Tag','MOVIE_REPEAT');
   set(h,'Enable',request);
end;

%--------------------------------------------

function gbT=StackSlider(gbT);

val = get(gcbo,'Value');
minVal = get(gcbo,'Min');
iVal = round(val);
gbT.movie.stackIndxVal = iVal-minVal+1;

indx = gbT.movie.stackIndxVal;

% update list
listH = findobj(gcbf,'Tag','STACK_LIST');
set(listH,'Value',gbT.movie.stackIndxVal);

% show the corresponding image
if(~isempty(gbT.viewPanelH))
   gbT.viewPanelH = ...
      uiViewPanelShowImg(gbT.data(:,:,indx),...
      0,gbT.viewPanelH);
else
   gbT.viewPanelH = ...
      uiViewPanelShowImg(gbT.data(:,:,indx),0);
end;

%--------------------------------------------

function gbT=StackList(gbT)
gbT.movie.stackIndxVal = get(gcbo,'Value');
indx = gbT.movie.stackIndxVal;

% update slider
sliderH = findobj(gcbf,'Tag','STACK_SLIDER');
if(strcmp(get(sliderH,'Enable'),'On'))
   set(sliderH,'Value',gbT.movie.stackIndx(indx));
end;

% show the corresponding image
if(~isempty(gbT.viewPanelH))
   gbT.viewPanelH = uiViewPanelShowImg(gbT.data(:,:,indx),...
      0,gbT.viewPanelH);
else
   gbT.viewPanelH = uiViewPanelShowImg(gbT.data(:,:,indx),0);
end;

%--------------------------------------------

function SetPath(gbT);
fileEditH = findobj(gcbf,'Tag','ACQ_FILE');
dirEditH = findobj(gcbf,'Tag','ACQ_DIR');
set(dirEditH,'String',gbT.acqDir);
set(fileEditH,'String',gbT.acqFile);

%--------------------------------------------

function EnableMovieCtls(gbT)
movieButtonH = findobj(gcbf,'Tag','MOVIE_BUTTON');
set(movieButtonH,'Enable','on');
if( isempty(gbT.movie.data))
	set(movieButtonH,'String',gbT.movie.buttonStr(1));
else
	set(movieButtonH,'String',gbT.movie.buttonStr(2));
end;   
repeatEditH = findobj(gcbf,'Tag','MOVIE_REPEAT');
set(repeatEditH,'String',int2str(gbT.movie.repeat));
set(repeatEditH,'Enable','on');
fpsEditH = findobj(gcbf,'Tag','FRAMES_SEC');
set(fpsEditH,'String',int2str(gbT.movie.fps));

%--------------------------------------------

function gbT=Browse(gbT)

initFileString = strcat(gbT.acqDir,filesep,'*.tif');

% get the new filename
[newFileString,newDirString] = ...
   uiputfile(initFileString,'Acquisition Panel: Set Filename');

if( isa(newFileString,'char') & isa(newDirString,'char') )
   % set the new data	
   gbT.acqDir = newDirString(1:length(newDirString)-1);
   gbT.acqFile = newFileString;
end;

%--------------------------------------------

function gbT=OpenFg(gbT)
% check if a framegrabber is active - close if yes
eval(['mexFgInterface(' 39 'close' 39 ')'],'');
gbT.fg.active=get(gcbo,'Value');
if (gbT.fg.active~=1)
   %enable grab buttons
   request = 'on';
	% value-2 for c++ starts counting at 0
	eval(['mexFgInterface(' 39 'open' 39 ',gbT.fg.active-2)'],['request=' 39 'off' 39 ';']);
else
   %diseable grab buttons
   request = 'off';
end;
acqFramSecH = findobj(gcbf,'Tag','FRAMES_SEC');
acqTimeH = findobj(gcbf,'Tag','ACQ_TIME');
shutTimeH = findobj(gcbf,'Tag','SHUTTER_TIME');
grabToViewH = findobj(gcbf,'Tag','GRAB_TO_VIEW');
grabToFileH = findobj(gcbf,'Tag','GRAB_TO_FILE');
grabStackH = findobj(gcbf,'Tag','GRAB_STACK');
brosweH = findobj(gcbf,'Tag','BROWSE');
fileEditH = findobj(gcbf,'Tag','ACQ_FILE');
dirEditH = findobj(gcbf,'Tag','ACQ_DIR');
set(acqFramSecH,'Enable',request);
set(acqTimeH,'Enable',request);
set(shutTimeH,'Enable',request);
set(grabToViewH,'Enable',request);
set(grabToFileH,'Enable',request);
set(grabStackH,'Enable',request);
set(brosweH,'Enable',request);
set(fileEditH,'Enable',request);
set(dirEditH,'Enable',request);

%--------------------------------------------

function SetShutterTime(gbT)
shutTimeH=findobj(gcbf,'Tag','SHUTTER_TIME');
sTime=str2num(get(shutTimeH,'String'));
if(isempty(gbT.subFrame))
   [s,dum]=mexfgInterface('grab',[1,1,2000,2000],1,sTime);
else
   [s,dum]=mexfgInterface('grab',gbT.subFrame,1,sTime);
end;


