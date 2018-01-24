function SetUserData (figHandle, object, replace, oName)
%SETUSERDATA connects any object to the figure 'figHandle'
%
% SYNOPSIS SetUserData(gHandle, object, oName)
%
% INPUT figHandle: a figure handle
%       object   : any object (number, string, struct...)
%       replace  : if replace = 1 the object is repalced in case of a name
%                  conflict. For any other value of replace an error occurs in
%                  case of a name conflict.
%       oName    : (optional) the name string under which the object is stored,
%                  by default the name of the 'object' variable is taken.
%
% c: 10/8/99	dT

uD = get(figHandle,'UserData');

% Check which name to use
if(nargin==4)
   name=oName;
else
   name=inputname(2);
end;

%Check if UserData is empty
if (isempty(uD))
   %initialize var
   uD=[];
else
   % Check for name conflict
   if ((any(strcmp(name,fieldnames(uD)))) & (replace~=1))
      error(['There exists already a field named ' 39 name 39]);
   end;
end;

uD=setfield(uD,name,object);
set(figHandle,'UserData',uD);
