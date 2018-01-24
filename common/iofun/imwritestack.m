function imwritestack(stack,firstfilename)
%IMWRITESTACK writes an image stack into a series of image files
%
% SYNOPSIS imwritestack(stack,firstfilename)
%
% INPUT stack      : image stack
%       firstfname : name of the first image to be written
%                    including the full path
%                    the actual filename must consist of 
%                    - alphanumeric body
%                    - numeric number
%                    - extension
%       
%
% OUTPUT none
%
% SEE ALSO imreadstack

if(ndims(stack) ~= 3)
   error('invalid stack entered');
end;
sDim = size(stack);

[fpath,fname,fno,fext]=getfilenamebody(firstfilename);

if(isempty(fname) | isempty(fno) | isempty(fext) )
   error('invalid first filename specified');
end;

if(isempty(fpath))
   fpath = pwd;
end;

nImg = sDim(3);

for(i=1:nImg)
   imIndx = str2num(fno)+i-1;
   imName= strcat(fpath,filesep,fname,int2str(imIndx),fext);
   imwrite(stack(:,:,i),imName,'Compression','none');
end;