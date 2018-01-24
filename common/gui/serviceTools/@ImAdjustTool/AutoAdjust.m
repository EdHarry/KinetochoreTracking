function AutoAdjust(imT)
% Adjusts gamma by setting the minimal value to 0 and the maximal to 1

cbFig = gcbf;

% Looking for stored image data
fig=findobj('Tag','UIVIEWPANEL');
figure(fig);
img=get(findobj(fig,'Type','image'),'CData');

% Setting min and max value of the image
mn=min(min(img));
mx=max(max(img));

% Using min and max to correct display 
newImT=ImAdjustTool(imT.figure);
newImT.newCMap=[];
newImT.oldInt=[mn mx];
newImT.newInt=[0 1];

% SetValues writes the new values into the text fileds IN1, IN2, OUT1, OUT2 (and GAMMA)
SetValues(newImT);
% SetContrast
SetContrast(newImT,cbFig);

% Functions

function SetValues(imT) 
cbFig = gcbf;
set(findobj(cbFig,'Tag','IN1'),'String',num2str(imT.oldInt(1)));
set(findobj(cbFig,'Tag','IN2'),'String',num2str(imT.oldInt(2)));
set(findobj(cbFig,'Tag','OUT1'),'String',num2str(imT.newInt(1)));
set(findobj(cbFig,'Tag','OUT2'),'String',num2str(imT.newInt(2)));
set(findobj(cbFig,'Tag','GAMMAEDIT'),'String',num2str(imT.gamma));
set(findobj(cbFig,'Tag','GAMMA'),'Value',min((imT.gamma)^.25,2));
 
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
fplot(['(' out2 '-' out1 ')*((x-' in1 ')/('in2 '-' in1 '))^' gam '+'  out1],[imT.oldInt imT.newInt]);
set(gca,'XGrid','On','YGrid','On','XLim',[0 1],'YLim',[0 1]);
hold on
plot([0 imT.oldInt(1)],[imT.newInt(1) imT.newInt(1)]);
plot([imT.oldInt(2) 1],[imT.newInt(2) imT.newInt(2)]);
hold off

% reset string
resetH = findobj(cbFig,'Tag','RESET');
set(resetH,'String','Reset');

% Display
displ('Intensity display adjusted. The actual values have not been modified.');
