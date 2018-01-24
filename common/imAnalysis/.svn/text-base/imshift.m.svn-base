function J = imshift(I,d,i0,md)
%IMSHIFT shifts an image by a vector
%
% The coordinate system is LEFT-HANDED with (1,1) being the center
% of the upper left corner pixel.
%
% SYNOPSIS J = imshift(I,d,i0,md)
%
% INPUT    I  : intensity image
%          d  : displacement vector [d(1), d(2)]
%               d(1) x1- (horizontal) shift
%               d(2) x2- (vertical) shift
%               the displacement is sepcified in a LEFT HANDED system 
%
%          i0 : dummy intensity for undefined values
%          md : method for interpolation (optional)
%               (cf. interp2)
% OUTPUT   J : image of size equal to I (optional)
%          undefined values are filled with 0 (dummy values)
%	        if no output argument is given the function will display
%          the map in the current axis    
%
% SEE ALSO imrotate, imresize, interp2
dim = size(I);

rem = d - fix(d);

if(any(rem ~= 0))
    [x1,x2]=meshgrid(1-d(1):1:dim(2)-d(1),1-d(2):1:dim(1)-d(2));
    
    if( nargin == 4 )
        J = interp2(I,x1,x2,md);
    else
        J = interp2(I,x1,x2);
    end;
    
    % replace NaN values by i0
    % mask = isnan(J);
    % J = setval(J,mask,i0);
    J(find(isnan(J))) = i0;
else
    % the displacement is an integer, thus just copy the matrix to the 
    % right place
    J = I;
    
    if d(2) <= 0
        J = [I(1-d(2):end,:); i0*ones(-d(2),dim(2))];
    elseif d(2) > 0	
        J = [i0*ones(d(2),dim(2)); I(1:end-d(2),:)];
    end
    
    if d(1) > 0
        J = [i0*ones(dim(1),d(1)), J(:,1:end-d(1))];
    elseif d(1) <= 0
        J = [J(:,-d(1)+1:end), i0*ones(dim(1),-d(1))];
end
    
    
end;
    
% display image if requested
if ( nargout == 0 )
   imshow(J);
end;
