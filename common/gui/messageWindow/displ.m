function displMsg(msg);

% Opens/activates display
msgH=findobj('Tag','DISPLAY_MSG');
if isempty(msgH)
   msgH=displayMsg;
end

% Looks for message string handle
txtH=findobj('Tag','MESSAGEFIELD');

% Displays string in message window
txtH=findobj('Tag','MESSAGEFIELD');
set(txtH,'String',msg);

