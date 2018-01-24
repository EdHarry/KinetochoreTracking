function gbT=SelectInteractive(gbT)

oldZoomStatus = 0;
oldPixelStatus = 0;

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
if(isempty(gbT.movie.stackIndx))
   % crop single image
   [gbT.data,gbT.subFrame] = imcrop;
   uiViewPanelShowImg(gbT.data,0,gbT.viewPanelH);
else
   % crop full stack
   [dummyI,gbT.subFrame] = imcrop;
   for( i = 1:length(gbT.movie.stackIndx) ) 
      dummyI(:,:,i) = imcrop(gbT.data(:,:,i),gbT.subFrame);
   end;
   listH = findobj(gcbf,'Tag','UIACQSTACKLIST');
   indx = get(listH,'Value');
   gbT.data=dummyI;
   gbT.viewPanelH = ...
      uiViewPanelShowImg(gbT.data(:,:,gbT.movie.stackIndxVal),...
      0,gbT.viewPanelH);
end;


%Set new roi to framegrabber (if an fg is active)- ignore errors
if (gbT.fg.active~=1)
   eval(['mexFgInterface(' 39 'grab' 39 ',gbT.subFrame)'],'');
end;

% set the old text back
set(textH,'String',oldText);
set(textH,'Visible',oldTextStatus);

% restore the old view panel status
if(oldPixelStatus)
   % toggle the pixel menu status back to on
   uiViewMenuPixelFnc(pixMenuH,figH);
else
   if(oldZoomStatus)
      % toggle the zoom menu status back to on
      uiViewMenuZoomFnc(zoomMenuH,figH);
   end;
end
