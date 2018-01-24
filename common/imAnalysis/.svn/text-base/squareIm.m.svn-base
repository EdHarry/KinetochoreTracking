function I = squareIm(fg,bg,side,ctr,dim,nse)
%SQUAREIM generation of a square image
% binary image with I = fg on background bg)
%
% SYNOPSIS I = squareIm(fg,bg,side,ctr,dim,nse)
%
% INPUT    fg : foreground intensity
%          bg : background intensity
%          side : side length
%          ctr: coordinates of the center
%               ctr(1) x- (horizontal) center coordinate
%               ctr(2) y- (vertical) center coordinate
%          dim: image dimension [nRows, nCols]
%               CAUTION: dim is compatible with the size() function
%                        but is NOT expressed in the left-handed 
%                        system.
%          nse: variance of the added zero mean Gaussian noise
% OUTPUT   I : (double) intensity image (optional)
%	        if no output argument is given the function will display
%          the map in the current axis

% check the compatibility of the parameters
sideEven = mod(side,2);
if(sideEven)
   error('side parameter NOT even');
end;
sideH = side/2;
if ( (ctr(1) - sideH < 1) | ...
      (ctr(1) + sideH > dim(2)) | ...
      (ctr(2) - sideH < 1) | ...
      (ctr(2) + sideH > dim(1)))
   error('imcompatible settings for ctr, side, dim');
end;

% create uniform bg image
I = bg * ones(dim);

% add fg rows
I((ctr(2)-sideH):(ctr(2)+sideH), ...
   (ctr(1)-sideH):(ctr(1)+sideH)) = fg;

% add noise
nse = abs(nse);
if(nse>0)
   I = imnoise(I,'gaussian',0,nse);
end;

% show the image if requested
if ( nargout == 0 )
   imshow(I);
end;


