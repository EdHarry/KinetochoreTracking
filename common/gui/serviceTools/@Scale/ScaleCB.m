function ScaleCB(sc)
%SCALE\SCALECB scale callback
%
% c: 12/8/99	dT

cbObjectH = gcbo;
cbFig =gcbf;

whoCalled = get(cbObjectH,'Tag');
switch whoCalled
	case 'SET'
      % Remove an exisiting scale from the figure
      if(any(ishandle(sc.lineH)))
         delete(sc.lineH);
         delete(sc.textH);
      end;
      % make scale figure active
      figure(sc.figure);
      % Get values
      origin = ginput(1);
      teH=findobj(cbFig,'Tag','POPCOLOR');
      sc.color = get(teH,'Value');
      switch(sc.color)
       case 1
         color=[1 1 1];
       case 2
         color=[0 0 0];
      end;
      teH=findobj(cbFig,'Tag','LENGTH');
      sc.length = str2num(get(teH,'String'));
      teH=findobj(cbFig,'Tag','PIXELSIZE');
		sc.pixelSize =str2num(get(teH,'String'));      
      teH=findobj(cbFig,'Tag','ADDTEXT');
      sc.text =get(teH,'String');
      len = sc.length/sc.pixelSize;
      sc.lineH = line([origin(1)-len/2 , origin(1)+len/2],[origin(2),origin(2)],...
         				'Color',color);
      sc.textH(1) = text(origin(1),origin(2)-7,sc.text,...
         				'Color',color,'HorizontalAlignment','center');   
     	% Connect scale to figure
   	SetUserData(sc.figure,sc,1,'myScale');
     	% Connect the Scale to the options window   
      set(gcbf,'UserData',sc);      
   case 'CLOSE'
      close(cbFig);
   case 'LENGTH'
      teH=findobj(cbFig,'Tag','ADDTEXT');
   	set(teH,'String',[get(gcbo,'String') ' \mum']);
end;