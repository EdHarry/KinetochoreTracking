function imwritestacknd(stack,firstfilename,depth)
%IMWRITESTACKND writes an image stack into a series of either 
% 8 or 16 bit image files
%
% SYNOPSIS imwritestack(stack,firstfilename,depth)
%
% INPUT stack      : image stack created by IMREADSTACKND
%                    (must be of class 'double' and normalized 
%                    to 0..1)
%       firstfname : name of the first image to be written
%                    including the full path
%                    the actual filename must consist of 
%                    - alphanumeric body
%                    - numeric number
%                    - extension
%       depth      : bit depth (8 or 16 bit)
%       
% OUTPUT none
%
% SEE ALSO IMREADSTACKND

% Check for input stack
if class(stack)=='uint8' | class(stack)=='uint16'
    imwritestack(stack,firstfilename);
    return;
end

if ~(class(stack)=='double' & max(max(max(stack)))<=1)
   error('The input stack is not a valid IMWRITESTACKND input parameter');
end

if(ndims(stack) ~= 3)
   error('invalid stack entered');
end;
sDim = size(stack);

if ~(depth==8 | depth==16)
   error('parameter ''depth'' must be either 8 or 16');
end

% Load
[fpath,fname,fno,fext]=getfilenamebody(firstfilename);

if(isempty(fname) | isempty(fno) | isempty(fext) )
   error('invalid first filename specified');
end;

if(isempty(fpath))
   fpath = pwd;
end;

nImg = sDim(3);

% Converting to 8 or 16 bit
switch depth
case 8
   stack=(2^8-1)*stack;
   stack=round(stack);
   stack=uint8(stack);
case 16
   stack=(2^16-1)*stack;
   stack=round(stack);
   stack=uint16(stack);
otherwise
end

% Writing image files
for(i=1:nImg)
   imIndx = str2num(fno)+i-1;
   imName= strcat(fpath,filesep,fname,int2str(imIndx),fext);
   imwrite(stack(:,:,i),imName,'Compression','none');
end;