function rgb = imOverlayOriMap(img,ori,mask,l,c)
%IMOVERLAYORIMAP overlays an orientation map to an image
%
% SYNOPSIS rgb = imOverlayOriMap(img,ori,mask,c)
%
% INPUT img : grey value image
%       ori : orientation ma [-pi,pi]
%       mask: binary image or image, where values == 1 control the drawing of an vector
%       l   : length of vector
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

if(max(max(mask))>1)
   error('invalid mask entered');
end;

auxMask = zeros(size(mask));
for( i = l:size(mask,1)-l )
   for( j = l:size(mask,2) - l)
      if(mask(i,j) == 1)
         xS = [j + l/2*cos(ori(i,j)),i+l/2*sin(ori(i,j))];
         xE = [j - l/2*cos(ori(i,j)),i-l/2*sin(ori(i,j))];
         xL = bresenham(xS,xE);
         for(k=1:size(xL,2))
            auxMask(xL(2,k),xL(1,k)) = 1;
         end;
      end;
   end;
end;

rgb = imOverlayMask(img,auxMask,c);



         
