function I = hanning(fg,bg,r0,ctr,dim,nse,ell)
%HANNING generation of a Hanning window 
% (spot with I = fg on background bg)
%
% The coordinate system is LEFT-HANDED with (1,1) being the center
% of the upper left corner pixel.
%
% SYNOPSIS I = hanning(fg,bg,r0,ctr,dim,nse,ell)
%
% INPUT    fg : foreground intensity
%          bg : background intensity
%          r0 : radius of hanning window 
%               (in case of elliptic window r0 = a )
%          ctr: coordinates of the center
%               ctr(1) x- (horizontal) center coordinate
%               ctr(2) y- (vertical) center coordinate
%          dim: image dimension [nRows, nCols]
%               CAUTION: dim is compatible with the size() function
%                        but is NOT expressed in the left-handed 
%                        system.
%          nse: variance of the added zero mean Gaussian noise
%          ell: row vector with elliptic parameters (optional)
%               ell = [ ecc, azi [rad] ], 
%               where ecc = 1 - b^2/a^2 ; 0 < ecc < 1;
%                     azi positive (clockwise) angle between 
%                         long ellipse axis and x- (horizontal) axis
% OUTPUT   I : (double) intensity image (optional)
%	        if no output argument is given the function will display
%          the map in the current axis    
c = (bg - fg)/2;
for row = 1:dim(1)
   for col = 1:dim(2)
      r = norm([row,col] - [ctr(2),ctr(1)]);
      if ( r > r0 )
         I(row,col) = bg;
      else
         I(row,col) = bg - c * ( 1 + cos (pi * r/r0) ); 
      end;
   end;
end;

% consider ellipse parameters if given
if( nargin > 6 )
   if ( (ell(1)>0) & (ell(1)<1))
      lambda = sqrt(1-ell(1)^2);
      A1 = [1,0 ; 0, lambda];
      A2 = [cos(ell(2)),-sin(ell(2)) ; sin(ell(2)), cos(ell(2))];
      J = imaffine(I,A2*A1,ctr,bg);
      I = J;
   end;
end;

% add noise
nse = abs(nse);
if(nse>0)
   I = imnoise(I,'gaussian',0,nse);
end;

% show the image if requested
if ( nargout == 0 )
   imshow(I);
end;

   