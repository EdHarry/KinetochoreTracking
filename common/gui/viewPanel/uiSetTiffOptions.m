function uiSetTiffOptions

% Sets all menu entries unchecked
par=get(gcbo,'Parent');
ch=get(par,'Children');
for counter1=1:length(ch)
   set(ch(counter1),'Checked','Off');
end
% Activates current menu option
set(gcbo,'Checked','On');

% Sets current adpeth and algorithm values
adal=get(gcbo,'UserData');

% Writes adepth and algorith values into 'UserData' of menu Load
ld=findobj(gcbf,'Tag','UIVIEWMENU_FILE_OPEN');
set(ld,'UserData',adal);




