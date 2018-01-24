function openFg
%OPENFG opens the globally specified framegrabber board
%
% The type of framegrabber is read from the global variable fgType__
%
% SYNOPSIS openFg
%
% INPUT none
%
% OUTPUT none
%

% version: started ??  GD
%

global fgType__;

% set the open command
if(strcmp(fgType__,'DT'))
   openCmd = 'dtFgOpen';
else
   return;
end;

% open the board
eval(openCmd);
