function I2 = imcropnint(I,rect,md)
%IMCROPDD crops an image at non-integer position
%
% SYNOPSIS I2 = imcropnint(I,rect,md)
%
% INPUT    I   : image
%          rect: [min(1) min(2) width height]
%                these values are specified in spatial coordinates
%          md  : interpoaltion method (cf. interp2)
%
% OUTPUT   I2  : cropped image of size width+1 height+1,
%                if rect is not inside I then I2 = []

[x,y]=meshgrid(rect(1):1:rect(1)+rect(3),rect(2):1:rect(2)+rect(4));

if( nargin == 4 )
   J = interp2(I,x,y,md);
else
   J = interp2(I,x,y);
end;

% check whether there are NaN's
mask = isnan(J);
if(max(max(mask)) == 1)
   I2 = [];
else
   I2 = J;
end;
