function closeFg
%CLOSEFG closes the globally specified framegrabber board
%
% The type of framegrabber is read from the global variable fgType__
%
% SYNOPSIS closeFg
%
% INPUT none
%
% OUTPUT none
%

global fgType__;

% set the close command
if(strcmp(fgType__,'DT'))
   closeCmd = 'dtFgClose';
else
   return;
end;

% open the board
eval(closeCmd);
