function convertStackFrameNames()
%CONVERTSTACKFRAMENAMES converts the frame names in a stack directory from the native (GD) to a sortable name format
%
% convertStackFrameNames()
%
% INPUT none : filename of the  first frame is specified interactively
%              the filename must consist of 
%					- alphanumeric body
%              - numeric number
%              - extension
%
% OUTPUT none : image stack
%        
% NOTE   the routine changes all the filenames which belong to a stack
%        series in a directory to the better sortable nomenclature 
%        with leading zeros. 
%        Example: test8.tif -> test008.tif
% 
%        The definition of how many laeding zero are inserted is taken 
%        from the total umber of frames in the stack.
%    
%        The old filenames are OVERWRITTEN

% started: 8-June-1999 GD

% declare variables
oldDir = [];

% get the  first filename
[fName,dirName] = uigetfile('*.tif','convertStackFramenames ...');

if( isa(fName,'char') & isa(dirName,'char'))
   firstfilename = [dirName,fName];
else
   return;
end;

[fpath,fname,fno,fext]=getfilenamebody(firstfilename);

if(isempty(fname) | isempty(fno) | isempty(fext) )
   error('invalid first filename specified');
end;

if(~isempty(fpath))
	% change to stack directory
	oldDir = cd(fpath);
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

imIndx = str2num(fno);
searchName = [fname,int2str(imIndx),fext];
if(~isempty(fileList))
   % get last index of the frame list
   while(~isempty(strmatch(lower(searchName),fileList)))
      imIndx = imIndx + 1;
      searchName= [fname,int2str(imIndx),fext];
   end;
   
   % define how many digits are necessary to represent the names
   maxIndx = imIndx-1;
   nDigits = floor(log10(maxIndx))+1;
   
   % go over all files between the start and the end index
   for(imIndx = str2num(fno):maxIndx)
      searchName= [fname,int2str(imIndx),fext];
      
      % look up table for how many leading zeros are necessary
      indxDigits = floor(log10(imIndx))+1;
      nAuxZeros = nDigits - indxDigits;
      switch nAuxZeros
      case 0, auxZeros = '';
      case 1, auxZeros = '0';
      case 2, auxZeros = '00';
      case 3, auxZeros = '000';
      case 4, auxZeros = '0000';
      case 5, auxZeros = '00000';
      otherwise, error('too large frame index');
      end;
      
      % move the filename (which is in MATLAB copy and delete)
      if(auxZeros > 0)
         newName = [fname,auxZeros,int2str(imIndx),fext];
         [status,msg]=copyfile(searchName,newName);
         if(status == 0)
            error(msg);
         else
            disp(['moved ',searchName,' to ',newName]);
            delete(searchName);
         end;
      end;
   end;
else 
   return
end;
   
% change back to original directory
if(~isempty(oldDir))
   cd(oldDir);
end;

