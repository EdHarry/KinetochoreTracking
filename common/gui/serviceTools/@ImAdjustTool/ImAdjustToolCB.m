function ImAdjustToolCB(imT)
%ImAdjustTool\ImAdjustToolCB  callback
%
% c: 17/8/99	dT

cbObjectH = gcbo;
cbFig =gcbf;

whoCalled = get(cbObjectH,'Tag');
switch whoCalled
	case{'IN1','IN2','OUT1','OUT2'}
      imT = GetInterval(imT);
      imT=SetContrast(imT,cbFig);   	
   case 'GAMMA'
      imT = GetInterval(imT);
      imT.gamma =get(gcbo,'Value')^4;
      gam = num2str(imT.gamma);
      %set new gamma value in tool window
      gamEdH=findobj(cbFig,'Tag','GAMMAEDIT');
      set(gamEdH,'String',gam);
      imT = SetContrast(imT,cbFig);
      
   case 'GAMMAEDIT'   
      imT = GetInterval(imT);
      gam = get(gcbo,'String');
      imT.gamma=str2num(gam);
      %set new gamma value in tool window
      gamSlidH=findobj(cbFig,'Tag','GAMMA');
     	set(gamSlidH,'Value',min((imT.gamma)^.25,2));
      imT = SetContrast(imT,cbFig);
      
   case 'CLOSE'
      close(cbFig);
      
   case 'RESET'
      if(strcmp(get(cbObjectH,'String'),'Last'))
         imT = SetContrast(imT,cbFig);
         set(cbObjectH,'String','Reset');
         SetValues(imT);
      else
         newImT=ImAdjustTool(imT.figure);
			newImT.newCMap=[];
			newImT.oldInt=[0 1];
			newImT.newInt=[0 1];
			newImT.gamma=1;
         SetValues(newImT);
        	SetContrast(newImT,cbFig);
			resetH = findobj(cbFig,'Tag','RESET');
         set(resetH,'String','Last');
      	% make image figure active
      	figure(imT.figure);      
         colormap(imT.oldCMap);
		end;
end;
% Connect the Scale to the options window
if(ishandle(imT.panelH))
	set(imT.panelH,'UserData',imT);
end;

function SetValues(imT)
cbFig = gcbf;   
set(findobj(cbFig,'Tag','IN1'),'String',num2str(imT.oldInt(1)));
set(findobj(cbFig,'Tag','IN2'),'String',num2str(imT.oldInt(2)));
set(findobj(cbFig,'Tag','OUT1'),'String',num2str(imT.newInt(1)));
set(findobj(cbFig,'Tag','OUT2'),'String',num2str(imT.newInt(2)));
set(findobj(cbFig,'Tag','GAMMAEDIT'),'String',num2str(imT.gamma));
set(findobj(cbFig,'Tag','GAMMA'),'Value',min((imT.gamma)^.25,2));

function imT = GetInterval(imT)
cbFig = gcbf;   
% Get values from imAdjustTool window
in1=get(findobj(cbFig,'Tag','IN1'),'String');
in2=get(findobj(cbFig,'Tag','IN2'),'String');
imT.oldInt=[str2num(in1) str2num(in2)];
out1=get(findobj(cbFig,'Tag','OUT1'),'String');
out2=get(findobj(cbFig,'Tag','OUT2'),'String');
imT.newInt=[str2num(out1) str2num(out2)];
   
function imT = SetContrast(imT,cbFig)
% make image figure active
figure(imT.figure);
%Back to strings
gam = num2str(imT.gamma);
in1= num2str(imT.oldInt(1));
in2= num2str(imT.oldInt(2));
out1= num2str(imT.newInt(1));
out2= num2str(imT.newInt(2));

imT.newCMap= imadjust(imT.oldCMap,imT.oldInt,imT.newInt,imT.gamma);
colormap(imT.newCMap);
      
% plot gamma curve: F(x)=(out2-out1)*((x-in1)/(in2-in1))^gamma+out2
figure(cbFig);
fplot(['(' out2 '-' out1 ')*((x-' in1 ')/(' in2 '-' in1 '))^' gam '+'  out1],[imT.oldInt imT.newInt]);
set(gca,'XGrid','On','YGrid','On','XLim',[0 1],'YLim',[0 1]);
hold on
plot([0 imT.oldInt(1)],[imT.newInt(1) imT.newInt(1)]);
plot([imT.oldInt(2) 1],[imT.newInt(2) imT.newInt(2)]);
hold off

% reset string
resetH = findobj(cbFig,'Tag','RESET');
set(resetH,'String','Reset');

  