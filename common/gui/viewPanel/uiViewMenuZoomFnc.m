function uiViewMenuZoomFnc(exH,exFigH);
% callback function
ownH = gcbo;
ownFigH = gcbf;

if(nargin > 0)
   ownH = exH;
   if(nargin == 2)
      ownFigH = exFigH;
   end;
end;

ownStatus = get(ownH,'Checked');
textH = findobj(ownFigH,'Type','uicontrol','Tag','UIVIEWTEXT');

if(strcmp(ownStatus,'off'))   	
   axesH = get(ownFigH,'CurrentAxes');
   imgH = findobj(axesH,'Type','image','Parent',axesH);
   data = get(imgH,'CData');
   if(~isempty(data))
	   pixMenuH = findobj(get(ownH,'Parent'),...
	   	'Tag','UIVIEWMENU_TOOLS_PIXEL');
      if(strcmp(get(pixMenuH,'Checked'),'on'))
         % toggle the pixel menu status back to off
         uiViewMenuPixelFnc(pixMenuH);
      end;
      set(ownH,'Checked','on');
      set(textH,'String','left button: zoom in; right button: zoom out');
      set(textH,'Visible','on');
     	zoom on;
   end;
else
   set(ownH,'Checked','off');
   set(textH,'Visible','off');
   zoom(ownFigH,'off');
end;
