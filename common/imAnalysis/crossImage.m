function I = crossImage(fg,bg,cross,ctr,dim,nse)
%CROSSIMAGE generation of a cross image
% binary image with I = fg on background bg)
%
% SYNOPSIS I = cross(fg,bg,cross,ctr,dim,nse)
%
% INPUT    fg : foreground intensity
%          bg : background intensity
%          cross : [ length, thickn ]
%                  both components must be even
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

% check the compatibility of the length, thickn, dim, and ctr
crossEven = mod(cross,2);
if(crossEven(1) | crossEven(2))
   error('cross parameters are NOT even');
end;
lengthH = cross(1)/2;
thicknH = cross(2)/2;
if ( (ctr(1) - lengthH < 1) | ...
      (ctr(1) + lengthH > dim(2)) | ...
      (ctr(2) - lengthH < 1) | ...
      (ctr(2) + lengthH > dim(1)))
   error('imcompatible settings for ctr, cross, dim');
end;

% create uniform bg image
I = bg * ones(dim);

% add fg colums
I((ctr(2)-lengthH):(ctr(2)+lengthH), ...
   (ctr(1)-thicknH):(ctr(1)+thicknH)) = fg;

% add fg rows
I((ctr(2)-thicknH):(ctr(2)+thicknH), ...
   (ctr(1)-lengthH):(ctr(1)+lengthH)) = fg;

% add noise
nse = abs(nse);
if(nse>0)
   I = imnoise(I,'gaussian',0,nse);
end;

% show the image if requested
if ( nargout == 0 )
   imshow(I);
end;


