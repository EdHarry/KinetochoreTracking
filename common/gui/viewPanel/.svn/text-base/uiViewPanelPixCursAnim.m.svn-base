function uiViewPanelPixCursAnim(action)
% service function to handle the cursor animation under 
% the tools_pixel menu
global y;
switch(action)
case 'start',
   axesH = get(gcbf,'CurrentAxes');
   imgH = findobj(axesH,'Type','image','Parent',axesH);
   currPt= get(axesH,'CurrentPoint');
   y(1)=round(currPt(1,1));
   y(2)=round(currPt(1,2));
   set(gcbf,'WindowButtonMotionFcn','uiViewPanelPixCursAnim move');
   set(gcbf,'WindowButtonUpFcn','uiViewPanelPixCursAnim stop');
   uiViewPanelPixCursAnim move;
case 'move',
   axesH = get(gcbf,'CurrentAxes');
   imgH = findobj(axesH,'Type','image','Parent',axesH);
   currPt= get(axesH,'CurrentPoint');
   imgData = get(imgH,'CData');
   x(1)= round(currPt(1,1));
   x(2)= round(currPt(1,2));
   dist(1)=x(1)-y(1);
   dist(2)=x(2)-y(2);
   p = impixel(imgData,x(1),x(2));
   if(sum(isnan(p)) == 0)
      if(isa(imgData,'uint8'))
   	  tString = sprintf('x(1) = %4d, x(2) = %4d, xd = %4d, yd = %4d, (%6d, %6d, %6d) ', ...
      		   x(1), x(2), dist(1), dist(2), p(1), p(2), p(3));
	   else	
   	   tString = sprintf('x(1) = %4d, x(2) = %4d, xd = %4d, yd = %4d, (%6.2f, %6.2f, %6.2f) ', ...
            x(1), x(2), dist(1), dist(2), p(1), p(2), p(3));
      end;
   	textH = findobj(gcbf,'Tag','UIVIEWTEXT');
      set(textH,'String',tString);
   end;
case 'stop'
   set(gcbf,'WindowButtonMotionFcn','');
   set(gcbf,'WindowButtonUpFcn','');
end;

   