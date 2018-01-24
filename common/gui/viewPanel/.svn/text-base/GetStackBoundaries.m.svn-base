function [gmin,gmax,fileMin,fileMax]=getStackBoundaries(firstfilename)
%GetStackBoundaries looks for global min and max values through a stack
%
% SYNOPSIS [gmin,gmax,fileMin,fileMax]=getStackBoundaries(firstfilename)
%
% INPUT firstfilename    : name of the first greyvalue image to be read
%                          including the full path
%                          the actual filename must consist of 
%                          - alphanumeric body
%                          - numeric number
%                          - extension
%
% OUTPUT gmin,gmax       : global min and max intensity values over the
%                          entire stack
%        fileMin, FileMax: files where global min and max have been found, respectively.
%
%
% NOTE  the search for filenames is case insensitive, i.e. the names in a stack
%       series can have different cases.


oldDir = [];

% Check for first file name validity
[fpath,fname,fno,fext]=getfilenamebody(firstfilename);

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

% Gets the list of all files and subdirs
dirListing = dir;

% get all relevant filenames (excludes dir names)
iEntry = 1;
fileList = {};
for( i = 1:length(dirListing))
   if(~dirListing(i).isdir)
      fileList(iEntry) = lower({dirListing(i).name});
      iEntry = iEntry + 1;
   end;
end;

% Only files with the given extension must be considered
nEntries = 0;
imIndx = str2num(fno);
l_fno=length(num2str(fno));

% First valid file name
searchName= [fname,num2str(imIndx,['%.' num2str(l_fno) 'd']),fext];

% Search for global min and max
if(~isempty(fileList))
   counter=0;
   while( ~isempty(strmatch(lower(searchName),fileList)))
      counter=counter+1;
      nEntries = nEntries + 1;
      index(nEntries) = imIndx;
      temp= imread(searchName);
      imIndx = imIndx + 1;
      % Looks for global min and max
      xmin=min(min(temp));
      xmax=max(max(temp));
      if counter==1
         gmin=xmin;
         gmax=xmax;
         fileMin=searchName;
         fileMax=searchName;
      else
         if gmin>=xmin
            gmin=xmin;
            fileMin=searchName;
         end
         if gmax<=xmax
            gmax=xmax;
            fileMax=searchName;
         end
      end
      clear temp;
		% Next file name
      searchName= [fname,num2str(imIndx,['%.' num2str(l_fno) 'd']),fext];
   end;
end;

% type conversion
gmin=double(gmin);
gmax=double(gmax);

% change back to original directory
if(~isempty(oldDir))
   cd(oldDir);
end;

