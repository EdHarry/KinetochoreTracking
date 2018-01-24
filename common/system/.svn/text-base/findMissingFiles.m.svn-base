function [missingList]=findMissingFiles(completeList,excludeDirs)
%findMissingFiles checks for essential Matlab files which are not found in the defined paths
%
%SYNOPSIS missingList=findMissingFiles(completeList)
%
%INPUT    completeList: cell array of functions on which program 'programName' depends, generated with the
%           command 'list=depfun(programName);' on a machine where the program works
%           If completeList is empty, a dialog box allows you to sear
%         excludeDirs : directories to exclude in the path search. The
%           match is case-sensitive, and the drive name is lower case.
%           If you have specified excludeDirs in the subfunction, you can
%           instead supply a numeric selection.
%
%OUTPUT   missingList: cell array of all function names that are not found in the matlab paths
%         
% note: with "for i=1:size(missingList,1);copyfile(missingList{i},'#folder#');end", you
%       can easily copy all missing files to the folder #folder#
%
%c: 1/03, Jonas Dorn
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%test input
if isempty(completeList)
    bait = uipickfiles('REFilter','\.m$');
    completeList = depfun(bait{:});
end
if ~iscell(completeList)
    error('sorry, wrong input. Generate cell array ''completeList'' with depfun');
end
if nargin < 2 || isempty(excludeDirs)
    excludeDirs = [];
elseif isnumeric(excludeDirs)
    excludeDirs = defaultExcludes(excludeDirs);
elseif ~iscell(excludeDirs)
    excludeDirs = {excludeDirs};
end
nExclude = numel(excludeDirs);


%initialize variables
k=1;
missingList{1,1}='';

%loop through completeList to build missingList
for i=1:size(completeList,1)
    %get first string
    currentFile=completeList{i};
    %extract filename
    fileSeparators=findstr(currentFile,filesep);
    %fileName: all chars after the last file separator. If no file separator, it is no file
    % also, check that the file is not on the matlab path
    if isempty(strmatch(matlabroot,currentFile)) &&  ~isempty(fileSeparators) 
        fileName=currentFile(fileSeparators(end)+1:end);
        
        %search paths for fileName
        isItThere=which(fileName, '-all');
        
        % check whether to exlcude dirs
        iter = 1;
        while ~isempty(isItThere) && iter <= nExclude
            % match all or beginning
            idx = strmatch(excludeDirs{iter},isItThere);
            if ~isempty(idx)
                isItThere(idx) = [];
            end
            iter = iter + 1;
        end
        
        %write fileName into output if not there
        if isempty(isItThere)
            missingList{k,1}=currentFile;
            k=k+1;           
        end
    end
end
missingList=sort(missingList);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% subfunctions
function excludeDirs = defaultExcludes(selection)
% defaultExcludes is a place to store default exclude directory lists

switch selection
    case 1
        excludeDirs = {'c:\data\jonas\matlab\common-tsri';...
            'c:\data\jonas\matlab\extern-tsri';...
            'c:\data\jonas\matlab\newFunctions';...
            'c:\data\jonas\matlab\mdxMisc';...
            'c:\data\jonas\matlab\chromdyn-tsri'};
    case 2
        excludeDirs = {'/Users/jonas/matlab'};
    otherwise
        excludeDirs = [];
end

