function B=nrm(A,d)
% nrm normalize a matrix to the range [0..1] (double), [0..255] (uint8) or [0..65535] (uint16)
%
% SYNOPSIS B=nrm(A,d)
%
% INPUT    A : matrix (image)
%          d : (optional: default = 1) 
%               bit depth (1, 8 or 16)
%                  1 - [0..1]     - output class: double
%                  8 - [0..255]   - output class: uint8
%                 16 - [0..65535] - output class: uint16
%
% OUTPUT   B : normalized A (with class change)

% By default, return a [0..1] double matrix
if nargin==1
    d=1;
end

% Normalize image to [0..1]
switch class(A)
case 'double'
   A=(A-min(A(:)))/(max(A(:))-min(A(:)));
otherwise
   A=double(A);
   A=(A-min(A(:)))/(max(A(:))-min(A(:)));
end   
% If A is complex (e.g. a transform)
B=abs(A);
% Conversion
switch d
case 1
case 8
   B=uint8((2^d-1)*B);
case 16
   B=uint16((2^d-1)*B);
otherwise
   error('The entered argument d is not valid');
end




