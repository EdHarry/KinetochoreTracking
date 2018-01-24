function img = imReadGui(type,title)
%IMREADGUI reads an image with a GUI / interface to imread()
%
% SYNOPSIS img = imReadGui(type)
%
% INPUT type: (optional) 'double' or 'uint8'; default -> 'double'
%             if set to 'struct' image becomes a structure with
%             *.data image matrix (double)
%             *.perm = 'M' permutation status set to MATLAB
%       title (optional) Title for the file dialogue. 
%              Default: 'imReadGui ...' 
%
% OUTPUT img : image matrix

if(nargin<1) || isempty(type)
   type = 'double';
end;
if nargin < 2 || isempty(title)
    title = 'imReadGui ...';
end

[fName,dirName] = uigetfile('*.tif',title);

if( isa(fName,'char') & isa(dirName,'char'))
   aux = imread([dirName,fName]);
else
   return;
end;

if(~strcmp(type,'uint8'))
   aux = double(aux)/255;
else
   img = aux;
   return;
end;

if(strcmp(type,'struct'))
   img.data = aux;
   img.perm = 'M';
else
   img = aux;
end;
