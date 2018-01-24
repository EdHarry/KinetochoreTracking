function uiViewMenuProfileFnc
% callback function

figH = gcbf;
textH = findobj(figH,'Type','uicontrol','Tag','UIVIEWTEXT');
set(textH,'String','left-click to select, right-click to finish');
set(textH,'Visible','on');
ImageProfile(figH);
set(textH,'Visible','off');
