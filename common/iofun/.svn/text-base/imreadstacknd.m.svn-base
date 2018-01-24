function [stack,msg,res,index,nEntries,outSubFrm,origColor]=imReadStackND(firstfilename,gmin,gmax,maxImg,inSubFrm)
% SYNOPSIS   [nstack,msg,res,index,nEntries,outSubFrm,origColor]=imReadStackND(firstfilename,gmin,gmax,maxImg,inSubFrm)
%
% INPUT  firstfilename: name of the first greyvalue image to be read
%                       including the full path
%                       the actual filename must consist of 
%                          - alphanumeric body
%                          - numeric number
%                          - extension
%        gmin, gmax   : intensity boundaries of the stack as returned by getStackBoundaries
%        maxImg       : (optional) maximal number of images to be processed
%        inSubFrm :     (optional) definition of subframe to be read
%                       [ul(1), ul(2), width, height]
%
% OUTPUT nstack       : normalized stack
%        msg,res      : return info on success or error (for panel)
%        index        : index number of images, e.g., if the first filename
%                       is 'test8.tif' the first index entry will be 8
%        nEntries     : number of filenames in this directory belonging to stack;
%                       for instance: directory contains 'test1.tif' - 'test100.tif'
%                                 maxImg is set to 50
%                           thus, length(index) == 50; nEntries == 100.
%        outSubFrm    : subframe definition which has actually been cropped
%                       [ul(1), ul(2), width, height]
%        origColor:     (optional) if this argument is supplied then the stack is kept in
%                       original color format (origColor = 'rgb' or 'gray').
%                       if it is omitted the stack is always converted to 'gray' format.


stack = [];
index = [];
oldDir = [];

% Check

if(nargin < 4)
   maxImg = inf;
end;

% firstfilename, gmin, gmax must be passed to the function
if (isempty(firstfilename) | isempty(gmin) | isempty(gmax))
   error('firstfilename, gmin, gmax must be passed to the function!');
end

% Load
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
   	if(nargin == 5)
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
   
   %allocate memory for stack which makes the loop go much faster
   height=outSubFrm(3);
   width=outSubFrm(4);
   while( ~isempty(strmatch(lower(searchName),fileList)) & (nEntries < maxImg))
      nEntries = nEntries + 1;
   end
%    if(nEntries == maxImg)
%       while(~isempty(strmatch(lower(searchName),fileList)))
%          nEntries = nEntries + 1;
%          imIndx = imIndx + 1;
%          searchName= [fname,num2str(imIndx,['%.' num2str(l_fno) 'd']),fext];
%       end;
%    end;

   if(nargout>=5) 
      if(isImgRGB(auxI))   %color images
         stack=zeros(width,height,3,nEntries);
      else           %monochrom image
         stack=zeros(width,height,nEntries);
      end
   end
   %nEntries = 0;
   
   %while( ~isempty(strmatch(lower(searchName),fileList)) & (nEntries < maxImg))
   for i_img=1:nEntries
      %nEntries = nEntries + 1;
      index(i_img) = imIndx;
      auxI = imread(searchName);
      % handle color images
      if(nargout>=5)
         if(isImgRGB(auxI))
            stack(:,:,:,i_img) = imcrop(auxI,outSubFrm);
			else
            stack(:,:,i_img) = imcrop(auxI,outSubFrm);
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
         stack(:,:,i_img) = imcrop(auxI,outSubFrm);
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

% Conversion of I to type double (necessary for normalization)
stack=double(stack);

% Normalization
stack=(stack-gmin)/(gmax-gmin);

% Check for normalization
ch=max(stack(:));
if ch>1
   msg='ERROR: Normalization failed! Select higher gmax and/or lower gmin';
   res=-1;
else
   msg=strcat('Normalization successful! [gmin and gmax are valid]');
   res=0;
end
   
% change back to original directory
if(~isempty(oldDir))
   cd(oldDir);
end;


