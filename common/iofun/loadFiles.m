function loadFiles(includeString,excludeString,directory,includeSubDirectories,selectionMode)
%LOADFILES is a utility to load a list of files. mat-files are loaded into the workspace, figures are displayed, and the rest is opened in the editor
%
%loadFiles is heavily based on searchFiles - hence the identical syntax
%
%SYNOPSIS loadFiles(includeString,excludeString,directory,includeSubDirectories,selectionMode)
%
%INPUT    includeString: string contained in the filenames you are looking for
%         excludeString (opt): string not contained in the filenames you are looking for
%         directory (opt): directory to search. if zero or 'pwd', current directory is searched
%                               if 'ask' (default), program asks for directory
%         includeSubDirectories (opt): whether to search subdirectories or not (0/{1})
%         selectionMode (opt): which file(s) to select if there are several files matching the
%                               includeString within one directory
%                              {'all'}-all; 'new'-newest; 'old'-oldest; 'GUI'-open GUI to select one file
%
%
%c: 11-03 jonas
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%---test input---
if nargin<1|isempty(includeString)
    error('not enough input arguments or empty includeString')
end

%includeString
if ~isstr(includeString)
    error('includeString has to be a string')
end

%excludeString
if nargin>1&~isempty(excludeString)
    if ~isstr(excludeString)
        error('excludeString has to be a string')
    end
else
    excludeString = [];
end

%directory (ask if necessary)
if nargin>2
    if isempty(directory)
        directory = 'ask';
    end
    switch directory
        case {0,'pwd'}
            directory = pwd;
        case {'ask'}
            directory = 'ask';
        otherwise
            if ~isdir(directory)
                error([directory,' is not a valid directory!'])
            end
    end
else
    directory = 'ask';
end

%includeSubDirectories
if nargin>3&~isempty(includeSubDirectories)
    if includeSubDirectories ~= 1
        includeSubDirectories = 0;
    end
else
    includeSubDirectories = 1;
end
if nargin<4|isempty(includeSubDirectories)
    includeSubDirectories = 1;
end

%selectionMode
if nargin<5|isempty(selectionMode)
    selectionMode = 'all';
else
    if ~(strcmp(selectionMode,'new')+strcmp(selectionMode,'old')+strcmp(selectionMode,'GUI')+strcmp(selectionMode,'all'))
        error('wrong selectionMode')
    end
end

%---end test input---


listOfFiles = searchFiles(includeString,excludeString,directory,includeSubDirectories,selectionMode);


%---load files
ct = 0;

for i = 1:size(listOfFiles,1)
    fileName = listOfFiles{i,1};
    %find fileExtension
    ptIdx = strfind(fileName,'.');
    if isempty(ptIdx)
        ptIdx = 0;
    end
    fileExt = fileName(ptIdx(end)+1:end);
    switch fileExt %currently recognizes mat-files, fig-files and m-files
        case 'mat'
            loadStruct = load([listOfFiles{i,2},filesep,fileName]);
            structFN = fieldnames(loadStruct);
            %write the data to the workspace
            for j = 1:length(structFN)
                varName = structFN{j};
                %check that this var does not exist yet
                if evalin('base',['exist(''',varName,''',''var'')'])
                    ct = ct+1;
                    varName = [varName,'_',num2str(ct)];
                end
                structFNjValue = eval(['loadStruct.',structFN{j},';']);
                assignin('base',varName,structFNjValue);
            end
        case 'fig'
            openfig([listOfFiles{i,2},filesep,fileName]);
        case 'm'
            edit([listOfFiles{i,2},filesep,fileName(1:end-2)]);
        otherwise %try to open in the editor
            try
            edit([listOfFiles{i,2},filesep,fileName]);
        catch
            disp(['could not open ',[listOfFiles{i,2},filesep,fileName],' with the editor.'])
        end
    end
           
end
%---end load files