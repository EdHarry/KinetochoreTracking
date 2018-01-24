function mHandle=myMessageBox(mHandle,message,name,offset)
% MYMESSAGEBOX display a message screen centered
% 
%
% SYNOPSIS mHandle=myMessageBox(handle,message,name,offset)
%
% INPUT mHandle : msgxboxHandle (if empty or zero a new box is opened)
%       message : string or string list
%       name    : msgbox name (only required on first call)
%       offset  : (optional) msgbox offset from screen center (in pixels [x y])
%
% OUTPUT mHandle : msgxboxHandle

% c: 7/3/03	dT

if nargin<4
    offset=[0 0];
end;

screenSize=get(0,'ScreenSize');
if isempty(mHandle) || mHandle==0  || ~ishandle(mHandle)
    if nargin<3
        name='';
    end;
    msbxSize=[100 100];
    mHandle=dialog('Position',[(screenSize(3)-msbxSize(1))/2+offset(1) (screenSize(4)+msbxSize(2))/2+offset(2) msbxSize(1) msbxSize(2)],'HandleVisibility','on','Name',name,'WindowStyle','normal','Tag','MessageBox');
    hText = uicontrol(mHandle,'Style', 'Text', 'String',message,'Tag','MESSAGE_TEXT');
else
    figure(mHandle);
    hText=findobj(mHandle,'Tag','MESSAGE_TEXT');
    set(hText,'String',message);
end;
msbxPos=get(mHandle,'Position');
extent=get(hText,'Extent');
msbxSize=extent(3:4)+20;
set(mHandle,'Position',[msbxPos(1:2) msbxSize(1) msbxSize(2)]);
% get actual msbSize
msbxPos=get(mHandle,'Position');
textPos=round([(msbxPos(3:4)-extent(3:4))/2 extent(3:4)]);
set(hText,'Position',textPos);
