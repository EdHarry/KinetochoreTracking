function gbT=GrabToStack(gbT)
% callback for the Grab Stack button in the acquisition panel

showCheckH = findobj(gcbf,'Tag','SHOW_ACQ');

% get timespan and frequency
timeEditH = findobj(gcbf,'Tag','ACQ_TIME');
fpsEditH = findobj(gcbf,'Tag','FRAMES_SEC');
gbT.movie.time = str2num(get(timeEditH,'String'));
gbT.movie.fps = str2num(get(fpsEditH,'String'));
if(isempty(gbT.movie.time) | (gbT.movie.time == 0))
   errordlg('no or invalid time specified','Acquisition Panel: Error');
   return;
end;
if(isempty(gbT.movie.fps) | (gbT.movie.fps == 0))
   errordlg('no or invalid frequency specified','Acquisition Panel: Error');
   return;
end;
nImg = ceil(gbT.movie.time*gbT.movie.fps) + 1;

if(nImg > gbT.movie.maxImgInStack)
   msg = sprintf('%d stack images requested; acquiring only %d',...
      nImg,gbT.movie.maxImgInStack);
   warndlg(msg,'Acquisition Panel: Error');
   nImg = gbT.movie.maxImgInStack;
   gbT.movie.time = (nImg-1)/gbT.movie.fps;
   set(timeEditH,'String',num2str(gbT.movie.time,'%8.1f'));
end;

% grab the stack
if(isempty(gbT.subFrame))
   gbT.data = grabSeq(nImg,gbT.movie.fps);
else
   gbT.data = grabSeq(nImg,gbT.movie.fps,gbT.subFrame);
end;

gbT.movie.stackIndx = 1:nImg;
gbT.movie.data = [];
gbT.movie.cMap = [];

% set parameters for UIACQSTACKSLIDER and UIACQSTACKLIST
sliderH = findobj(gcbf,'Tag','STACK_SLIDER');
listh = findobj(gcbf,'Tag','STACK_LIST');

set(listh,'String',...
   int2str(reshape(gbT.movie.stackIndx,nImg,1)), ...
   'Value',1);

set(sliderH,'Min',gbT.movie.stackIndx(1),...
   'Max',gbT.movie.stackIndx(nImg),...
   'SliderStep',[1/(nImg-1),1],...
   'Value',gbT.movie.stackIndx(1));

if(get(showCheckH,'Value'))
   if(~isempty(gbT.viewPanelH))
      gbT.viewPanelH = ...
         uiViewPanelShowImg(gbT.data(:,:,1),0,gbT.viewPanelH);
   else
      gbT.viewPanelH = ...
         uiViewPanelShowImg(gbT.data(:,:,1),0);
   end;
end;