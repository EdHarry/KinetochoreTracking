function rgb = imOverlayMask(img,mask,c)
%IMOVERLAYMASK overlays a mask on a grayvalue image and produces an rgb map
%
% SYNOPSIS rgb = imOverlayMask(img,mask,c)
%
% INPUT img : grey value image
%       mask: binary image or image, where values == 1 do mask the original image
%       c   : (optional) color code; default 'r'
%
% OUTPUT rgb : rgb map
%
% SEE ALSO colorCode2rgb for color code definitions

if(nargin < 3)
   c = 'r';
end

if(isa(img,'uint8'))
   img = double(img);
end;

if(isa(mask,'uint8'))
   mask = double(mask);
end;

if(max(max(mask))>1)
   error('invalid mask entered');
end;

maskI = ~mask;
aux = img .* maskI;
colVec = colorCode2rgb(c);
rgb = cat(3,aux,aux,aux)+cat(3,mask.*colVec(1),mask.*colVec(2),mask.*colVec(3));
