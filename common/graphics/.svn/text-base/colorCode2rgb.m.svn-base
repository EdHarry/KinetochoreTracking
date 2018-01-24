function rgbVec = colorCode2rgb(c)
%COLORCODE2RGB converts a color code to an rgb vector
%
% The following colors are supported:
%  'y' yellow
%  'm' magenta 
%  'c' cyan
%  'r' red
%  'g' green
%  'b' blue
%  'w' white
%  'k' black
%
% SYNOPSIS rgbVec = colorCode2rgb(c)
%
% INPUT c : color code
% 
% OUTPUT rgbVec : vector with the rgb value

switch c
case 'y', rgbVec = [1,1,0];
case 'm', rgbVec = [1,0,1];
case 'c', rgbVec = [0,1,1];
case 'r', rgbVec = [1,0,0];
case 'g', rgbVec = [0,1,0];
case 'b', rgbVec = [0,0,1];
case 'w', rgbVec = [1,1,1];
case 'k', rgbVec = [0,0,0];
otherwise, rgbVec = [], error('unknown color code')
end;
