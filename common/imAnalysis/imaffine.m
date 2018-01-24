function J = imaffine(I,A,ctr,i0,md)
%IMAFFINE warps the image according to affine trasformation.
% 
% The coordinate system is LEFT-HANDED with (1,1) being the center
% of the upper left corner pixel.
%
% By setting an appropriate matrix A and transformation center ctr
% also congruent and rigid body transformations can be performed.
% The transformation does not include coordinate center shifts. To
% achieve this use the function imshift()
%
% SYNOPSIS J = imaffine(I,A,ctr,i0,md)
%
% INPUT    I  : image
%          A  : 2x2 transformation matrix
%          ctr: transformation center
%               ctr(1) x- (horizontal) center coordinate
%               ctr(2) y- (vertical) center coordinate
%          i0 : dummy intensity for undefined values
%          md : method for interpolation (optional)
%               (cf. interpol2)
% OUTPUT   J : image of size equal to I (optional)
%	        if no output argument is given the function will display
%          the map in the current axis    
%
% SEE ALSO imrotate, imresize, imshift, interpol2

tUint8 = 0;

if isa(I,'uint8')
    I = double(I)/255;
    i0 = double(i0)/255;
    tUint8 = 1;
end;

dimI = size(I);
dimA = size(A);

% prepare inverse transformation
if( (dimA(1) ~= 2) | (dimA(2) ~=2) )
   error('Transformation matrix A must be 2x2 invertible');
end;
invA = inv(A);

[x,y] = meshgrid(1:dimI(2),1:dimI(1));


% transform the coordinates 
XI = invA *([(x(:)-ctr(1))';(y(:)-ctr(2))']);
xI = reshape(XI(1,:),dimI(1),dimI(2))+ctr(1);
yI = reshape(XI(2,:),dimI(1),dimI(2))+ctr(2);

% inverse transformation to prepare grid coordinates
%for(row = 1:dimI(1))
%   for(col = 1:dimI(2))
%      pos = [col;row] - [ctr(1); ctr(2)];
%      x(row,col) = invA(1,:) * pos + ctr(1);
%      y(row,col) = invA(2,:) * pos + ctr(2);
%   end;
%end;

% interpolate
if( nargin == 5 )
   J = interp2(I,xI,yI,md);
else
   J = interp2(I,xI,yI);
end;

% replace NaN values by i0
J(find(isnan(J))) = i0;


% transform class back to uint8 if necessary
if tUint8
    J = uint8(J*255);
end;

% display image if requested
if ( nargout == 0 )
   imshow(J);
end;
