function out = ReadSunras(filename)
% Read image from file in SUN raster format
%
% SYNOPSIS out = ReadSunras(filename)
%
% INPUT:  filename String
% OUTPUT: out uint8 array
%
% NOTICE: the function is a wrapper to the mexReadSunras function 
%         which must be placed in the private directory below this.

if ~(nargin == 1)
   error('1 input argument requested');
end;
if ~(nargout == 1)
   error('1 output argument requested');
end;
if ~ischar(filename)
   error('argument is not a string');
end;

out = mexReadSunras(filename);
