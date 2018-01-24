function obj = GetUserData (figHandle, oName, remove)
%GETUSERDATA gets any object connected to the figure 'figHandle'
%
% SYNOPSIS object = GetUserData(gHandle, oName, remove)
%
% INPUT  figHandle: a figure handle
%        oName    : the name string under which the object is stored
%        remove   : (optional) if remove = 1 then the connection of 'oName' 
%                   to the figure is removed
%
% OUTPUT object   : the object
%
% c: 12/8/99	dT

uD = get(figHandle,'UserData');

% Check if name already exists in the list
if (~isempty(uD) & (any(strcmp(oName,fieldnames(uD)))))
   obj = getfield(uD,oName);
   % Remove the field?
   if ((nargin==3) & (remove==1))
      uD = rmfield(uD,oName);
      set(figHandle,'UserData',uD);
   end;
else
   obj = [];
end;
