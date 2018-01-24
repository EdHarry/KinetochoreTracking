function uiViewMenuPixelFnc(exH,exFigH)
% callback function
ownH = gcbo;
ownFigH = gcbf;

if(nargin > 0)
   ownH = exH;
   if(nargin == 2)
      ownFigH = exFigH;
   end;
end;

status = get(ownH,'Checked');
textH = findobj(ownFigH,'Type','uicontrol','Tag','UIVIEWTEXT');
axesH = get(ownFigH,'CurrentAxes');
imgH = findobj(axesH,'Type','image','Parent',axesH);

if(strcmp(status,'off'))
   data = get(imgH,'CData');
   if(~isempty(data))
      zoomMenuH = findobj(get(ownH,'Parent'),...
	   	'Tag','UIVIEWMENU_TOOLS_ZOOM');
      if(strcmp(get(zoomMenuH,'Checked'),'on'))
         % toggle the zoom menu status back to off
         uiViewMenuZoomFnc(zoomMenuH);
      end;	
      set(ownH,'Checked','on');
      set(textH,'String','press left button');
      set(textH,'Visible','on');
      set(ownFigH,'Pointer','fullcrosshair');
      set(imgH,'ButtonDownFcn','uiViewPanelPixCursAnim start');
   end;
else
   set(ownH,'Checked','off');
   set(textH,'Visible','off');
   set(ownFigH,'Pointer','arrow');
   set(imgH,'ButtonDownFcn','');
end;
