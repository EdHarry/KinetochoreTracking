function rgb = mkColorImg(chan);
%MKCOLORIMG makes a color image from an intensity file
%
% SYNOPSIS rgb = mkColorImg(chan)
% 
% INPUT chan : channel to which the intensity image shold go 
%              either 'r', 'g', or 'b'
%
% OUTPUT rgb : rgb image 
%
% NOTE : if a color image is selected interactively, then the routine 
%        will convert it to an intensity image first, and use this as the
%        requested channel 
%

% STARTED 30-5-2000 GC

% get file interactively
i = imreadGui;

if isempty(i)
   error('invalid image read');
end;

dim = ndims(i);

switch dim
case 2, rgb = constructRgb(i,chan);
case 3,
   % make sure that a 3 channel image is an intensity image indeed
   rgb = constructRgb(rgb2gray(i),chan);
otherwise
   error('invalid image format');
end;


%-----------------------------------------------------------------


function rgb = constructRgb(grey,chan)
% constructor of the color image

r = zeros(size(grey));
g = zeros(size(grey));
b = zeros(size(grey));

switch chan
case 'r', r = grey;
case 'g', g = grey;
case 'b', b = grey;
otherwise 
   error('invalid channel specification');
end;

rgb = cat(3,r,g);
rgb = cat(3,rgb,b);
