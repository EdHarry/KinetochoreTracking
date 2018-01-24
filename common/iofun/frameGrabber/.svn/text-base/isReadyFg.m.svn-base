function ans = isReadyFg
%ISREADYFG checks whether the globally specified framegrabber board
%         is ready to be used
%
% The type of framegrabber is read from the global variable fgType__
%
% SYNOPSIS openFg
%
% INPUT none
%
% OUTPUT 1 if true
%
global fgType__;

% set the ready command
if(strcmp(fgType__,'DT'))
   readyCmd = 'dtFgIsOpen';
else
   ans = 0;
   return;
end;

% evaluate the command
ans = eval(readyCmd);
