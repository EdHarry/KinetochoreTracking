function imwritend(image,filename,depth)
%IMWRITEND writes an image to a graphics file.
%
%   SYNOPSIS IMWRITEND(filename,image,depth)
%
%   INPUT   image:      image (of class 'double') loaded and 
%                       normalized to [0..1] by IMREADND  
%           filename:   name of the file with entire path
%           depth:      bit depth (8 or 16 bit)
%       
% OUTPUT none
%
% SEE ALSO IMREADND

% *************************
%
% INPUT PARAMETER CHECK
%
% *************************

% Check for existance/validity of parameter depth
if isempty(depth)
   error('A value for ''depth'' must be defined');
end
%
if ~(depth==8 || depth==16)
   error('parameter ''depth'' must be either 8 or 16');
end

% Check for input image
if ~(isa(image,'double') && max(max(image))<=1)
   error('The input image is not a valid IMWRITEND input parameter');
end

% Setting image to either 8 or 16 bit
image=(2^depth-1)*image;
image=round(image);

% Defining class
switch depth
case 8
   image=uint8(image);
case 16
   image=uint16(image);
otherwise
end

% Writing file
imwrite(image,filename,'Compression','none');


