function [stack,index,nEntries,outSubFrm,origColor]=imreadstack(firstfilename,maxImg,inSubFrm)
%IMREADSTACK reads a series of image files (grey value images) into an image stack
%
% SYNOPSIS [stack,index,nEntries,outSubFrm]=imreadstack(firstfilename,maxImg,inSubFrm)
%
% INPUT firstfname : name of the first greyvalue image to be read
%                    including the full path
%                    the actual filename must consist of 
%							- alphanumeric body
%                    - numeric number
%                    - extension
%       maxImg : (optional) maximal number of images
%       inSubFrm : (optional) definition of subframe to be read
%                  [ul(1), ul(2), width, height]
%
% OUTPUT stack : image stack
%        index : index number of images, e.g., if the first filename
%                is 'test8.tif' the first index entry will be 8
%        nEntries : number of filenames in this directory belonging to stack;
%                   for instance: directory contains 'test1.tif' - 'test100.tif'
%                                 maxImg is set to 50
%                           thus, length(index) == 50; nEntries == 100.
%        outSubFrm : subframe definition which has actually been cropped
%                  [ul(1), ul(2), width, height]
%        origColor: (optional) if this argument is supplied then the stack is kept in
%                   original color format (origColor = 'rgb' or 'gray').
%                   if it is omitted the stack is always converted to 'gray' format.
%
% NOTE  the search for filenames is case insensitive, i.e. the names in a stack
%       series can have different cases.

stack = [];
index = [];
oldDir = [];

if(nargin == 1)
   maxImg = inf;
end;

[fpath,fname,fno,fext]=getFilenameBody(firstfilename);

if(isempty(fname) | isempty(fno) | isempty(fext) )
   error('invalid first filename specified');
end;

if(~isempty(fpath))
	% change to stack directory
   oldDir = cd(fpath);
else
   %check if it is in the matlab search path
   tempName=which(firstfilename);
   if(~isempty(tempName))
      [fpath,fname,fno,fext]=getfilenamebody(tempName);
      oldDir = cd(fpath);
	end;
end;

dirListing = dir;

% get all relevant filenames
iEntry = 1;
fileList = {};
for( i = 1:length(dirListing))
   if(~dirListing(i).isdir)
      fileList(iEntry) = lower({dirListing(i).name});
      iEntry = iEntry + 1;
   end;
end;

nEntries = 0;
imIndx = str2num(fno);
l_fno=length(num2str(fno));
searchName= [fname,num2str(imIndx,['%.' num2str(l_fno) 'd']),fext];

if(~isempty(fileList))
   % check the validity of subframe specification and reset it if necessary
   if(~isempty(strmatch(lower(searchName),fileList)))
      auxI = imread(searchName);
      [height,width,depth] = size(auxI);
      frm = [ 1, 1, width, height];
   	if(nargin == 3)
         outSubFrm = isinside(inSubFrm,frm,1);
      else
         outSubFrm = frm;
      end;
      %check color type
      if(nargout>=5)
      	if(isImgRGB(auxI))
            origColor = 'rgb';
         else
            origColor = 'gray';
         end;
      end;
   end;
   while( ~isempty(strmatch(lower(searchName),fileList)) & (nEntries < maxImg))
      nEntries = nEntries + 1;
      index(nEntries) = imIndx;
      auxI = imread(searchName);
      % handle color images
      if(nargout>=5)
         if(isImgRGB(auxI))
            stack(:,:,:,nEntries) = imcrop(auxI,outSubFrm);
			else
            stack(:,:,nEntries) = imcrop(auxI,outSubFrm);
         end;
      else
         if(isImgRGB(auxI))
            %Check if it is really a rgb or just a single channel rgb
            singleRGB=double(any(any(auxI(:,:,1)))) + double(any(any(auxI(:,:,2))))...
               + double(any(any(auxI(:,:,3))));
            if (singleRGB==1)
               auxI=bitor(auxI(:,:,1),bitor(auxI(:,:,2),auxI(:,:,3)));
            else   
               auxI=rgb2gray(auxI);
            end;           
         end;
         stack(:,:,nEntries) = imcrop(auxI,outSubFrm);
      end;
      imIndx = imIndx + 1;
      searchName= [fname,num2str(imIndx,['%.' num2str(l_fno) 'd']),fext];
   end;
   if(nEntries == maxImg)
      while(~isempty(strmatch(lower(searchName),fileList)))
         nEntries = nEntries + 1;
         imIndx = imIndx + 1;
         searchName= [fname,num2str(imIndx,['%.' num2str(l_fno) 'd']),fext];
      end;
   end;
end;

% change back to original directory
if(~isempty(oldDir))
   cd(oldDir);
end;

